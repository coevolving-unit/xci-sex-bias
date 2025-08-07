library(biomaRt)
library(rcompanion)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(patchwork)
library(egg)
library(RColorBrewer)
library(tibble)
library(data.table)
library(stringr)
library(mashr)

#### load and combine data ####

hsap = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
map = getBM(mart = hsap, attributes = c('ensembl_gene_id', 'external_gene_name','chromosome_name','start_position','end_position'))
saveRDS(map, file = 'map.rds')

map = readRDS('map.rds')

# ASE from Gylemo et al.

ase = read.csv('/data/gylemo.csv')
colnames(ase)[3] = 'external_gene_name'
colnames(ase)[6] = 'tissue'
ase = merge(ase, map, by = 'external_gene_name', all.x = T) 
ase$tissue = as.factor(ase$tissue)
levels(ase$tissue) = c('AdiposeSubcutaneous','AdiposeVisceralOmentum','AdrenalGland','ArteryAorta',
                       'ArteryCoronary','ArteryTibial','BrainAnteriorcingulatecortexBA24','BrainCaudatebasalganglia',
                       'BrainCerebellarHemisphere','BrainSpinalcordcervicalc1','Breast','ColonTransverse','EsophagusMucosa','EsophagusMuscularis',
                       'HeartAtrialAppendage','HeartLeftVentricle','KidneyCortex','Liver','Lung','MuscleSkeletal',
                       'Ovary','Pancreas','SkinSunExposedLowerleg','SkinNotSunExposedSuprapubic','Spleen','Stomach',
                       'Thyroid','Uteris','Vagina','Whole Blood')
ase %>% group_by(new_category) %>% summarise(unique_genes = n_distinct(external_gene_name))
ase$current_escape = ase$new_category

ase$current_escape = ifelse(ase$allelic_expression < 0.4, 'escape', 'inactive')
table(ase$new_category, ase$current_escape)

colnames(ase)
ase = subset(ase, chromosome_name == 'X')

check = ase %>% group_by(individual, tissue, external_gene_name) %>% summarise(n=n())
check = subset(check, n > 1)

write.csv(ase, file = 'ase.csv')

# sex bias from DeCasien et al.

sb = readRDS('/data/mashr_sex.rds')
beta = get_pm(sb)
beta = reshape2::melt(beta)
lfsr = get_lfsr(sb)
lfsr = reshape2::melt(lfsr)
sb = cbind(beta, lfsr$value)
colnames(sb) = c('ensembl_gene_id','tissue','beta','lfsr')
sb$bias = ifelse(sb$lfsr < 0.05 & sb$beta > 0, 'male', 'ns')
sb$bias = ifelse(sb$lfsr < 0.05 & sb$beta < 0, 'female', sb$bias)
table(sb$bias)
sb$tissue = as.factor(sb$tissue)
levels(sb$tissue)

sb = merge(sb, map, by = 'ensembl_gene_id', all.x = T) 
sb = subset(sb, chromosome_name %in% c(1:22,"X","Y"))

write.csv(sb, file = 'sb.csv')

# combine

co = merge(ase, sb, by = c('ensembl_gene_id','external_gene_name','tissue'))
co$tissue = droplevels(co$tissue)
levels(co$tissue)
table(co$tissue, co$new_category)
table(co$bias, co$new_category)
table(co$bias, co$current_escape)

co$cat3 = co$current_escape
co$cat3 = ifelse(co$new_category=='PAR','PAR',co$current_escape)
table(co$bias, co$cat3)

tissue.list = levels(co$tissue)
short.listnow = c('Adipose(Sub)','Adipose(Vis)','Adrenal','Artery(Aor)','Artery(Cor)','Artery(Tib)','Brain(BA24)','Brain(Caud)','Brain(Cblm)','Spinal(C1)','Colon(Trans)','Esoph(Muc)','Esoph(Musc)','Heart(Atr)','Heart(Ven)','Kidney(Cor)','Liver','Lung','Muscle','Pancreas','Skin(Sun)','Skin(NoSun)','Spleen','Stomach','Thyroid')
palette_Dark2 <- colorRampPalette(brewer.pal(17, "Dark2"))
tissue.colors = palette_Dark2(25)
tissue.shapes = c(0:24)

check = co %>% group_by(individual, tissue, external_gene_name) %>% summarise(n=n())
check = subset(check, n > 1)
View(check)

co2 <- co %>%
  group_by(external_gene_name, tissue) %>%
  summarize(
    current_escape = if (n_distinct(current_escape) > 1) "variable" else unique(current_escape),
    new_category = unique(new_category),
    bias = unique(bias),
    .groups = "drop"
  )

unique(subset(co2, new_category == 'PAR' & current_escape != 'escape')$external_gene_name)
View(subset(co, current_escape == 'inactive' & new_category == 'PAR'))
View(subset(co2, current_escape == 'inactive' & new_category == 'PAR'))
View(subset(co2, current_escape == 'variable' & new_category == 'PAR'))

write.csv(co, file = 'co.csv')
write.csv(co2, file = 'co2.csv')

#### end ####

#### Figure S2A ####

table(co2$current_escape, co2$new_category)
m = table(co2$current_escape, co2$new_category)
m
m/sum(m)*100
m/rowSums(m)*100
m / matrix(colSums(m), nrow = nrow(m), ncol = ncol(m), byrow = TRUE)*100
chisq.test(m)
pairwiseNominalIndependence(m, fisher = T, gtest = T, method = "bonferroni")
pairwiseNominalIndependence(t(m), fisher = T, gtest = T, method = "bonferroni")

# Create df_plot and clean up data as before
df_plot <- as.data.frame.matrix(m) %>%
  rownames_to_column("current_escape") %>%
  pivot_longer(-current_escape, names_to = "new_category", values_to = "count")

df_plot <- df_plot %>%
  mutate(new_category = dplyr::recode(new_category,
                               "escape_across_tissues" = "broadly escape",
                               "escape_data_for_single_tissue" = "escape 1 tissue",
                               "inactive_across_tissues" = "broadly inactive",
                               "inactive_data_for_single_tissue" = "inactive 1 tissue",
                               "PAR" = "PAR",
                               "variable_across_tissues" = "variable"))
df_plot$new_category <- factor(df_plot$new_category, 
                               levels = c("PAR", 
                                          "broadly escape", 
                                          "escape 1 tissue", 
                                          "variable", 
                                          "broadly inactive", 
                                          "inactive 1 tissue"))
current_escape_colors <- c(
  "escape" = "#56B4E9",     # Sky blue
  "inactive" = "#E69F00",    # Orange
  "variable" = "#F0E442"           # Yellow
)

# Compute counts and proportions
df_bar <- df_plot %>%
  group_by(current_escape, new_category) %>%
  summarise(n = sum(count), .groups = "drop") %>%
  group_by(new_category) %>%
  mutate(prop = n / sum(n),
         label_y = cumsum(prop) - (prop / 2))  # for label position

# Plot for counts (bar chart)
count_plot <- ggplot(df_bar, aes(x = new_category, y = n, fill = current_escape)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Count", x = "Gene-level category", fill = 'Tissue-level category') +
  scale_fill_manual(values = current_escape_colors) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x axis title for the count plot
        axis.text.x = element_blank(),   # Remove x axis labels for the count plot
        axis.ticks.x = element_blank(),  # Remove x axis ticks for the count plot
        legend.position = "none")        # Remove legend from count_plot

# Plot for proportions (stacked bar chart)
prop_plot <- ggplot(df_bar, aes(x = new_category, y = prop, fill = current_escape)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Proportion", x = 'Gene-level category', fill = 'Tissue-level category') +
  scale_fill_manual(values = current_escape_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +  # Apply to prop_plot only
  guides(fill = guide_legend(title = "Tissue-level category", nrow = 1))  # Ensure only one legend

# Combine both plots into one layout, sharing the x-axis and one legend
count_plot + prop_plot + plot_layout(nrow = 2, heights = c(1, 1))

#### end ####

#### Figure S2C #### 

table(subset(co2, new_category == 'PAR')$bias, subset(co2, new_category == 'PAR')$current_escape)
m = table(subset(co2, new_category == 'PAR')$bias, subset(co2, new_category == 'PAR')$current_escape)
m
m/sum(m)*100
m/rowSums(m)*100
m / matrix(colSums(m), nrow = nrow(m), ncol = ncol(m), byrow = TRUE)*100
chisq.test(m)
pairwiseNominalIndependence(m, fisher = T, gtest = T, method = "bonferroni")
pairwiseNominalIndependence(t(m), fisher = T, gtest = T, method = "bonferroni")

tbl <- matrix(c(197, 8,
                58, 3),
              nrow = 2, byrow = TRUE)
rownames(tbl) <- c("male", "others")
colnames(tbl) <- c("escape", "not_escape")
tbl
fisher.test(tbl)

df_plot <- as.data.frame.table(m) %>%
  rename(bias = Var1, current_escape = Var2, count = Freq)

df_bar <- df_plot %>%
  group_by(current_escape, bias) %>%
  summarise(n = sum(count), .groups = "drop") %>%
  group_by(current_escape) %>%
  mutate(prop = n / sum(n),  # Proportions within each current_escape category
         label_y = cumsum(prop) - (prop / 2)) 

bias_colors <- c(
  "female" = "#D55E00",   # Vermillion
  "male" = "#009E73",     # green blue
  "ns" = "#666"        # grey
)

count_plot <- ggplot(df_bar, aes(x = current_escape, y = n, fill = bias)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Count", x = "Tissue-level category", fill = 'Sex-bias') +
  scale_fill_manual(values = bias_colors) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x axis title for the count plot
        axis.text.x = element_blank(),   # Remove x axis labels for the count plot
        axis.ticks.x = element_blank(),  # Remove x axis ticks for the count plot
        legend.position = "none")        # Remove legend from count_plot

# Plot for proportions (stacked bar chart)
prop_plot <- ggplot(df_bar, aes(x = current_escape, y = prop, fill = bias)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Proportion", x = 'Tissue-level category', fill = 'Sex-bias') +
  scale_fill_manual(values = bias_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +  # Apply to prop_plot only
  guides(fill = guide_legend(title = "Sex-bias", nrow = 1))  # Ensure only one legend

# Combine both plots into one layout, sharing the x-axis and one legend
count_plot + prop_plot + plot_layout(nrow = 2, heights = c(1, 1))

#### end ####

#### Figure S2E #### 

table(subset(co2, new_category != 'PAR')$bias, subset(co2, new_category != 'PAR')$current_escape)
m = table(subset(co2, new_category != 'PAR')$bias, subset(co2, new_category != 'PAR')$current_escape)
m
m/sum(m)
m/rowSums(m)
m / matrix(colSums(m), nrow = nrow(m), ncol = ncol(m), byrow = TRUE)*100
chisq.test(m)
pairwiseNominalIndependence(m, fisher = T, gtest = T, method = "bonferroni")
pairwiseNominalIndependence(t(m), fisher = T, gtest = T, method = "bonferroni")

tbl <- matrix(c(240, 355,
                115, 2800),
              nrow = 2, byrow = TRUE)
rownames(tbl) <- c("female", "others")
colnames(tbl) <- c("escape", "not_escape")
tbl
fisher.test(tbl)

df_plot <- as.data.frame.table(m) %>%
  rename(bias = Var1, current_escape = Var2, count = Freq)

df_bar <- df_plot %>%
  group_by(current_escape, bias) %>%
  summarise(n = sum(count), .groups = "drop") %>%
  group_by(current_escape) %>%
  mutate(prop = n / sum(n),  # Proportions within each current_escape category
         label_y = cumsum(prop) - (prop / 2)) 

bias_colors <- c(
  "female" = "#D55E00",   # Vermillion
  "male" = "#009E73",     # green blue
  "ns" = "#666"        # grey
)

count_plot <- ggplot(df_bar, aes(x = current_escape, y = n, fill = bias)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Count", x = "Tissue-level category", fill = 'Sex-bias') +
  scale_fill_manual(values = bias_colors) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),  # Remove x axis title for the count plot
        axis.text.x = element_blank(),   # Remove x axis labels for the count plot
        axis.ticks.x = element_blank(),  # Remove x axis ticks for the count plot
        legend.position = "none")        # Remove legend from count_plot

# Plot for proportions (stacked bar chart)
prop_plot <- ggplot(df_bar, aes(x = current_escape, y = prop, fill = bias)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Proportion", x = 'Tissue-level category', fill = 'Sex-bias') +
  scale_fill_manual(values = bias_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +  # Apply to prop_plot only
  guides(fill = guide_legend(title = "Sex-bias", nrow = 1))  # Ensure only one legend

# Combine both plots into one layout, sharing the x-axis and one legend
count_plot + prop_plot + plot_layout(nrow = 2, heights = c(1, 1))

#### end ####

#### Figure S2B/D/E S4C ####

plot_data <- co2 %>%
  arrange(new_category, external_gene_name) %>%
  mutate(gene_factor = factor(external_gene_name, levels = unique(external_gene_name)))

ggplot(plot_data, aes(x = tissue, y = gene_factor, fill = current_escape)) +
  geom_tile(color = "white") +
  facet_grid(new_category ~ ., scales = "free", space = "free") +
  scale_fill_manual(values = current_escape_colors) +
  theme_article(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold")
  ) +
  labs(x = "Genes", fill = "Tissue-level category")

ggplot(plot_data, aes(x = tissue, y = gene_factor, fill = bias)) +
  geom_tile(color = "white") +
  facet_grid(new_category ~ ., scales = "free", space = "free") +
  scale_fill_manual(values = bias_colors) +
  theme_article(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 4),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0, face = "bold")
  ) +
  labs(x = "Genes", fill = "Tissue-level category")

# PAR genes + full escape + sex bias
pe = subset(co, PAR == 'PAR' & allelic_expression < 0.01)
ggplot(pe, aes(x = beta, y = reorder(external_gene_name, -beta), color = bias, shape = tissue)) +
  geom_point(size = 3) +
  scale_shape_manual(values = tissue.shapes) +
  scale_color_manual(values = bias_colors) +
  theme_article() +
  theme(axis.title.y = element_blank())

#### end ####

#### Figure S4A/B/D/E 1B S2F 1C ####

# plot all

co$mix = paste(co$bias, co$current_escape)
ggplot(co, aes(x = allelic_expression, y = beta, color = mix, shape = mix)) +
  geom_point(size = 3, alpha = 0.8) + 
  geom_smooth(method = 'lm', color = 'black') +
  geom_hline(yintercept = 0) +
  xlim(c(0,0.5)) +
  ylim(c(-1,1)) +
  theme_article()

# plot male-biased PAR + female-biased NPX

pl1 = subset(co, new_category != 'PAR' & bias == 'female')
pl2 = subset(co, new_category == 'PAR' & bias == 'male')
pl = rbind(pl1, pl2)

ggplot(pl, aes(x = allelic_expression, y = beta, color = current_escape, fill = bias, shape = bias)) +
  geom_point(aes(color = current_escape, fill = bias, shape = bias), size = 3, stroke = 1) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_shape_manual(values = c(21, 22)) +
  geom_vline(xintercept = 0.4, linetype = 'dashed') +
  scale_fill_manual(values = bias_colors) +
  scale_color_manual(values = current_escape_colors) +
  geom_smooth(aes(group = interaction(current_escape, bias)),method = 'lm',color = "black",fill = "grey70",alpha = 0.3) + 
  xlim(c(0,0.5)) +
  ylim(-1.1, 1.1) +
  theme_article()

ggplot(pl, aes(x = 1-allelic_expression, y = beta)) +
  geom_point(aes(color = tissue, shape = tissue), size = 3, stroke = 1) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_shape_manual(values = tissue.shapes) +
  geom_vline(xintercept = 0.4, linetype = 'dashed') +
  scale_color_manual(values = tissue.colors) +
  geom_smooth(aes(group = interaction(current_escape, bias)),method = 'lm',color = "black",fill = "grey70",alpha = 0.3) + 
  xlim(c(0.5,1)) +
  ylim(-1.1, 1.1) +
  theme_article() 

ggplot(pl, aes(y = beta, x = current_escape)) +
  facet_wrap(~bias) + 
  geom_boxplot(aes(color = current_escape, fill = bias)) + 
  scale_fill_manual(values = bias_colors) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ylim(-1.1, 1.1) + # excludes 1 outlier
  scale_color_manual(values = current_escape_colors) +
  theme_article()

s = subset(pl, current_escape == 'escape')
s$ae = 1-s$allelic_expression
mod = lm(beta ~ ae * PAR, data = s)
summary(mod)

# tests

m = matrix(table(co$bias, co$cat3), ncol=3, nrow=3)
rownames(m) = c('female','male','ns')
colnames(m) = c('escape','inactive','PAR')
m
m = m[-3,]
m
m/sum(m)
m/rowSums(m)
sweep(m, 2, colSums(m), "/")

chisq.test(m)
pairwiseNominalIndependence(m, fisher = T, gtest = T, method = "bonferroni")
pairwiseNominalIndependence(t(m), fisher = T, gtest = T, method = "bonferroni")

# within tissues

c = pl %>% group_by(tissue, current_escape, PAR) %>% summarise(n=n())
c = subset(c, current_escape == 'escape' & n > 3)
pl2 = subset(pl, tissue %in% c$tissue)

o = data.frame()
for(i in 1:length(unique(pl2$tissue))){
  np = subset(pl2, tissue == unique(pl2$tissue)[i] & current_escape == 'escape' & PAR == 'PAR')
  cp = cor.test(np$beta, 1-np$allelic_expression, method = 'pearson')
  nx = subset(pl2, tissue == unique(pl2$tissue)[i] & current_escape == 'escape' & PAR == 'nonPAR')
  cx = cor.test(nx$beta, 1-nx$allelic_expression, method = 'pearson')
  o[i,1] = np$tissue[1]
  o[i,2] = 'PAR'
  o[i,3] = cp$estimate
  o[i,4] = cp$p.value
  o[i,5] = dim(np)[1]
  o[i,6] = 'nonPAR'
  o[i,7] = cx$estimate
  o[i,8] = cx$p.value
  o[i,9] = dim(nx)[1]
}
View(o)
write.csv(o, file = 'tissue-specific-results.csv')

# per gene

c = pl %>% group_by(external_gene_name, current_escape) %>% summarise(n=n())
c = subset(c, current_escape == 'escape' & n > 3)
pl2 = subset(pl, external_gene_name %in% c$external_gene_name)

o = data.frame()
for(i in 1:length(unique(pl2$external_gene_name))){
  n = subset(pl2, external_gene_name == unique(pl2$external_gene_name)[i] & current_escape == 'escape')
  c = cor.test(n$beta, 1-n$allelic_expression, method = 'pearson')
  o[i,1] = n$external_gene_name[1]
  o[i,2] = n$PAR[1]
  o[i,3] = c$estimate
  o[i,4] = c$p.value
  o[i,5] = dim(n)[1]
}
View(o)
write.csv(o, file = 'gene-specific-results.csv')

pl2f = subset(pl2, bias == 'female' & external_gene_name %in% subset(o, V4 < 0.05)$V1)
ggplot(pl2f, aes(x = 1-allelic_expression, y = beta, color = bias, shape = tissue)) +
  geom_point(aes(color = bias, shape = tissue), size = 3, stroke = 1) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  scale_color_manual(values = bias_colors) +
  scale_shape_manual(values = 0:24) +
  facet_wrap(~external_gene_name) + 
  geom_smooth(aes(group = current_escape),method = 'lm',color = "black",fill = "grey70",alpha = 0.3) + 
  theme_article()+
  theme(legend.position = 'right')
pl2f = subset(pl2, PAR == 'nonPAR')
ggplot(pl2f, aes(x = 1-allelic_expression, y = beta, color = bias, shape = tissue)) +
  geom_point(aes(color = bias, shape = tissue), size = 3, stroke = 1) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  scale_color_manual(values = bias_colors) +
  scale_shape_manual(values = 0:24) +
  facet_wrap(~external_gene_name, scales = 'free') + 
  geom_smooth(aes(group = current_escape),method = 'lm',color = "black",fill = "grey70",alpha = 0.3) + 
  theme_article()+
  theme(legend.position = 'right')

pl2m = subset(pl2, bias == 'male' & external_gene_name %in% subset(o, V4 < 0.05)$V1)
ggplot(pl2m, aes(x = 1-allelic_expression, y = beta, color = bias, shape = tissue)) +
  geom_point(aes(color = bias, shape = tissue), size = 3, stroke = 1) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  ylim(c(0,1)) +
  scale_color_manual(values = bias_colors) +
  scale_shape_manual(values = 0:24) +
  facet_wrap(~external_gene_name) + 
  geom_smooth(aes(group = current_escape),method = 'lm',color = "black",fill = "grey70",alpha = 0.3) + 
  theme_article()+
  theme(legend.position = 'right')
pl2m = subset(pl2, PAR == 'PAR')
ggplot(pl2m, aes(x = 1-allelic_expression, y = beta, color = bias, shape = tissue)) +
  geom_point(aes(color = bias, shape = tissue), size = 3, stroke = 1) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  ylim(c(0,1)) +
  scale_color_manual(values = bias_colors) +
  scale_shape_manual(values = 0:24) +
  facet_wrap(~external_gene_name) + 
  geom_smooth(aes(group = current_escape),method = 'lm',color = "black",fill = "grey70",alpha = 0.3) + 
  theme_article()+
  theme(legend.position = 'right')

#### end ####

#### Figure 1A ####

p = subset(co, PAR == "PAR")[,c('ensembl_gene_id','allelic_expression','beta',
                                'current_escape','tissue','individual')]
k = merge(p, map, by = 'ensembl_gene_id', all.x = T)
k = subset(k, chromosome_name == 'X')
k = unique(k)
k = subset(k, external_gene_name != 'SPRY3')
k = subset(k, external_gene_name != 'XG')
k = k[order(k$start_position),]
k

tads = fread('/data/A-172_GSE147123_tad.bed', header = FALSE)
tads = subset(tads, V1 == 'chrX')

k$group = ifelse(k$start_position < tads$V3[1], 'A', 'X')
k$group = ifelse(k$start_position > tads$V3[2], 'B', k$group)
k$group = ifelse(k$start_position > tads$V3[3], 'C', k$group)
k$group = ifelse(k$start_position > tads$V3[4], 'D', k$group)
k$group = ifelse(k$start_position > tads$V3[5], 'E', k$group)
table(table(k$external_gene_name))
unique(k[,c('group','external_gene_name','start_position')])

unique_genes_kb <- k %>%
  mutate(start_kb = start_position / 1000)
x_min <- 10
x_max <- 2781.479

ggplot(unique_genes_kb) +
  geom_hline(yintercept = 0, linewidth = 2, color = "darkgreen") +
  geom_segment(aes(x = start_kb, xend = start_kb, y = 0, yend = 1), color = NA, linewidth = 0.5) +
  geom_text(aes(x = start_kb, y = 0.01, label = external_gene_name), 
            angle = 90, hjust = 0, size = 3) +
  geom_vline(xintercept = tads$V3[c(1:6)]/1000) +
  scale_x_continuous(breaks = seq(x_min, x_max, by = 100)) +
  labs(x = "Genomic Position (kp)", y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())

means <- k %>%
  group_by(group) %>%
  summarise(mean_ae = mean(allelic_expression, na.rm = TRUE))
midpoints <- k %>%
  group_by(group) %>%
  summarise(
    mean_ae = mean(allelic_expression, na.rm = TRUE),
    mid_gene = external_gene_name[which.min(abs(rank(start_position) - median(rank(start_position))))]
  )
midpoints <- midpoints %>%
  left_join(k %>% dplyr::select(group, external_gene_name, start_position), 
            by = c("group", "mid_gene" = "external_gene_name")) %>%
  distinct(group, .keep_all = TRUE)

ggplot(k, aes(x = reorder(external_gene_name, start_position), y = allelic_expression)) +
  ylim(c(0, 0.5)) +
  facet_grid(~group, scales = 'free_x', space = 'free_x') +
  geom_boxplot() + 
  geom_point(data = midpoints, 
             aes(x = reorder(mid_gene, start_position), y = mean_ae), 
             inherit.aes = FALSE, shape = 21, fill = "red", size = 3) + 
  theme_article()

mod = aov(1-allelic_expression ~ group, data = k)
summary(mod)
TukeyHSD(mod)

# Compute pairwise distances across all samples
k_pairs <- expand.grid(i = 1:nrow(k), j = 1:nrow(k)) %>%
  filter(i < j) %>%
  mutate(
    ae_i = k$allelic_expression[i],
    ae_j = k$allelic_expression[j],
    group_i = k$group[i],
    group_j = k$group[j],
    distance = abs(ae_i - ae_j),
    same_group = group_i == group_j
  )
k_pairs_named <- k_pairs %>%
  mutate(
    gene_i = k$external_gene_name[i],
    gene_j = k$external_gene_name[j]
  ) %>%
  filter(gene_i != gene_j) 

# Observed effect: difference in median distance
within <- k_pairs_named %>% filter(same_group) %>% pull(distance)
between <- k_pairs_named %>% filter(!same_group) %>% pull(distance)
median(between)
median(within)
mean(between)
mean(within)
observed_effect <- median(between) - median(within)

# Permutation test: shuffle group labels and recalculate
set.seed(123)
n_perm <- 1000
perm_effects <- replicate(n_perm, {
  shuffled_groups <- sample(k$group)
  k_pairs$group_i <- shuffled_groups[k_pairs$i]
  k_pairs$group_j <- shuffled_groups[k_pairs$j]
  k_pairs$same_group <- k_pairs$group_i == k_pairs$group_j
  
  w <- k_pairs %>% filter(same_group) %>% pull(distance)
  b <- k_pairs %>% filter(!same_group) %>% pull(distance)
  median(b) - median(w)
})

# P-value: how often permuted effects exceed observed
p_value <- mean(perm_effects >= observed_effect)

# Output
cat("Observed effect (median[between] - median[within]):", round(observed_effect, 4), "\n")
cat("Permutation test p-value:", p_value, "\n")

#### end ####

#### Figure S5 ####

pl1 = subset(co, new_category != 'PAR' & bias == 'female')
pl2 = subset(co, new_category == 'PAR' & bias == 'male')

k1 = pl1 %>% group_by(tissue) %>% summarise(ae = median(1-allelic_expression), beta = median(-beta))
k1$rankae = rank(k1$ae)
k1$rankb = rank(k1$beta)
ggplot(k1, aes(x = rankae, y = rankb)) +
  geom_point(size  =NA) +
  geom_smooth(method = 'lm') +
  geom_point(aes(x = rankae, y = rankb, color = tissue, shape = tissue)) +
  scale_shape_manual(values=c(0:24)) +
  theme_classic() +
  xlab('XCI escape rank') + ylab('female-biased expression rank')
cor.test(k1$rankae, k1$rankb, method = 'spearman')

k2 = pl2 %>% group_by(tissue) %>% summarise(ae = median(1-allelic_expression), beta = median(beta))
k2$rankae = rank(k2$ae)
k2$rankb = rank(k2$beta)
ggplot(k2, aes(x = rankae, y = rankb)) +
  geom_point(size  =NA) +
  geom_smooth(method = 'lm') +
  geom_point(aes(x = rankae, y = rankb, color = tissue, shape = tissue)) +
  scale_shape_manual(values=c(0:24)) +
  theme_classic() +
  xlab('XCI escape rank') + ylab('male-biased expression rank')
cor.test(k2$rankae, k2$rankb, method = 'spearman')

#### end ####

#### Figure S6 ####

his = read.csv('/data/hou-et-al-h3k27ac.csv')
colnames(his)[18] = 'external_gene_name'
his$tissue = as.factor(his$tissue)
levels(his$tissue) = c('BrainAnteriorcingulatecortexBA24','HeartAtrialAppendage','Lung','MuscleSkeletal')
write.csv(his, file = 'his.csv')

his2 = merge(co, his, by = c('tissue','external_gene_name'))
table(his2$sex.biase.direction, his2$bias, his2$current_escape, his2$PAR)
p = his2 %>% group_by(sex.biase.direction, bias, current_escape, PAR) %>% summarise(n=n())
p = subset(p, bias !='ns')

ggplot(p, aes(x = sex.biase.direction, y = bias, fill = n)) +
  geom_tile() + 
  geom_text(aes(label = n), color = "black", size = 3) +
  facet_grid(current_escape ~ PAR, space = 'free', scales = 'free') +
  theme_article() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = "Sex biased enhancer", y = "Sex biased expression", fill = "Count")

tbl <- matrix(c(41, 2, 14, 53), nrow = 2, byrow = TRUE,
              dimnames = list(GeneGroup = c("PAR_escape", "Other"),
                              H3K27acBias = c("MaleBiased", "FemaleBiased")))
tbl
sum(tbl)
fisher.test(tbl)

tbl <- matrix(c(24, 0, 44, 42), nrow = 2, byrow = TRUE,
              dimnames = list(
                GeneGroup = c("nonPAR_escape", "Other"),
                FemaleBiased = c("Yes", "No")
              ))
tbl
sum(tbl)
fisher.test(tbl)

#### end ####

#### Figure S1 #### 

# Oliva et al.

sb = read.csv('/data/oliva.csv')
colnames(sb) = c('ensembl_gene_id','external_gene_name','chr','beta','sd','lfsr','tissue')
sb$ensembl_gene_id = str_sub(sb$ensembl_gene_id,1,15)
sb$beta = -1*sb$beta
sb$bias = ifelse(sb$lfsr < 0.05 & sb$beta > 0, 'male', 'ns')
sb$bias = ifelse(sb$lfsr < 0.05 & sb$beta < 0, 'female', sb$bias)
table(sb$bias)
sb$tissue = as.factor(sb$tissue)
levels(sb$tissue) = c('AdiposeSubcutaneous','AdiposeVisceralOmentum','AdrenalGland','ArteryAorta',
                      'ArteryCoronary','ArteryTibial','Breast','BrainAnteriorcingulatecortexBA24','BrainCaudatebasalganglia',
                      'BrainCerebellarHemisphere','BrainSpinalcordcervicalc1','ColonSigmoid','ColonTransverse','EsophagusGastroesophagealJunction',
                      'EsophagusMucosa','EsophagusMuscularis',
                      'HeartAtrialAppendage','HeartLeftVentricle','KidneyCortex','Liver','Lung','MuscleSkeletal',
                      'Pancreas','Pituitary','SkinNotSunExposedSuprapubic','SkinSunExposedLowerleg','SmallIntestineTerminalIleum','Spleen','Stomach','Thyroid','Whole Blood')
sb = subset(sb, external_gene_name != 'XIST')

# DeCasien et al. 

sb2 = readRDS('/data/mashr_sex.rds')
beta = get_pm(sb2)
beta = reshape2::melt(beta)
lfsr = get_lfsr(sb2)
lfsr = reshape2::melt(lfsr)
sb2 = cbind(beta, lfsr$value)
colnames(sb2) = c('ensembl_gene_id','tissue','beta','lfsr')
sb2$bias = ifelse(sb2$lfsr < 0.05 & sb2$beta > 0, 'male', 'ns')
sb2$bias = ifelse(sb2$lfsr < 0.05 & sb2$beta < 0, 'female', sb2$bias)
table(sb2$bias)
sb2$tissue = as.factor(sb2$tissue)
levels(sb2$tissue)

# combine
co = merge(sb, sb2, by = c('ensembl_gene_id','tissue'))

# plot
ggplot(co, aes(x = beta.x, y = beta.y, color = tissue)) +
  geom_point(alpha = 0.8) +
  xlab('Oliva et al. 2020') +
  ylab('DeCasien et al. 2025') +
  geom_smooth(method = 'lm', alpha = 0.8) +
  theme_article()

t = unique(co$tissue)
o = data.frame()
for(i in 1:length(t)) {
  n = subset(co, tissue == t[i])
  c = cor.test(n$beta.x, n$beta.y, method = 'spearman')
  o[i,1] = t[i]
  o[i,2] = c$estimate
  o[i,3] = c$p.value
  o[i,4] = dim(n)[1]
}
o
write.csv(o,'oliva-vs-decasien.csv')

#### end ####

#### Figure S3 ####

sc = read.csv('/data/other-ASE-results.csv')
# Gylemo = abs(0.5 - (Xa/total))
# sclinax (AIDA, Japanese, XXY, tabula-sapiens, mutliome) = Xi:total
# SanRoman = Xi:Xa -> adjust
# Garieri = Xi:total
# Tukiainen = Xa:total -> adjust 

# averaged + per sample
# Garieri
# Tukiainen
# San Roman

# adjusted + raw
# San Roman

sc = sc[!(is.na(sc$xci.class)),]
length(unique(sc$external_gene_name))
table(sc$xci.class)
table(sc$dataset)
table(sc$celltype.organ.sample)

ggplot(sc, aes(x = ratio.xi.expression), alpha = 0.5) +
  facet_grid(dataset~xci.class) +
  geom_boxplot() +
  geom_vline(xintercept = 0.5, color = 'red') +
  theme_article()

sc <- sc %>%
  mutate(ratio.xi.expression = if_else(
    dataset %in% c("Tukanian_averaged","Tukiainen_persample") & !is.na(ratio.xi.expression),
    1-ratio.xi.expression, 
    ratio.xi.expression))

sc <- sc %>%
  mutate(ratio.xi.expression = if_else(
    dataset %in% c("SanRoman_averaged","SanRoman_persample") & !is.na(ratio.xi.expression),
    ratio.xi.expression / (1 + ratio.xi.expression),
    ratio.xi.expression))

ggplot(sc, aes(x = ratio.xi.expression), alpha = 0.5) +
  facet_grid(dataset~xci.class) +
  geom_vline(xintercept = 0.5, color = 'red') +
  geom_boxplot() +
  theme_article()

ggplot(sc, aes(x = reorder(external_gene_name, ratio.xi.expression), y = ratio.xi.expression)) +
  geom_boxplot() +
  facet_grid(dataset~xci.class, scales = 'free_x', space = 'free') +
  theme_article() +
  #ylim(c(0.2, 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3),
        axis.title.x = element_blank())

ggplot(sc, aes(x = ratio.xi.expression, y = celltype.organ.sample)) +
  geom_boxplot() +
  facet_grid(dataset~xci.class, scales = 'free_y', space = 'free') +
  theme_article() +
  #ylim(c(0.2, 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3),
        axis.title.x = element_blank())

ggplot(sc, aes(x = reorder(external_gene_name, ratio.xi.expression), y = ratio.xi.expression)) +
  geom_boxplot() +
  facet_grid(.~xci.class, scales = 'free_x', space = 'free') +
  theme_article() +
  #ylim(c(0.2, 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3),
        axis.title.x = element_blank())

scpar = subset(sc, xci.class == 'PAR1')
scpar = subset(scpar, measure != "raw")
table(scpar$celltype.organ.sample)
scpar = merge(scpar, subset(map, chromosome_name == 'X'), by = 'external_gene_name')
table(scpar$dataset)

#View(subset(scpar, ratio.xi.expression > 0.5))

ggplot(scpar, aes(x = reorder(external_gene_name, start_position), y = ratio.xi.expression)) +
  geom_boxplot() +
  theme_article() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

tads = fread('/data/A-172_GSE147123_tad.bed', header = FALSE)
tads = subset(tads, V1 == 'chrX')

scpar$group = ifelse(scpar$start_position < tads$V3[1], 'A', 'X')
scpar$group = ifelse(scpar$start_position > tads$V3[2], 'B', scpar$group)
scpar$group = ifelse(scpar$start_position > tads$V3[3], 'C', scpar$group)
scpar$group = ifelse(scpar$start_position > tads$V3[4], 'D', scpar$group)
scpar$group = ifelse(scpar$start_position > tads$V3[5], 'E', scpar$group)
table(scpar$external_gene_name, scpar$group)

means <- scpar %>%
  group_by(group) %>%
  summarise(mean_ae = mean(ratio.xi.expression, na.rm = TRUE))
midpoints <- scpar %>%
  group_by(group) %>%
  summarise(
    mean_ae = mean(ratio.xi.expression, na.rm = TRUE),
    mid_gene = external_gene_name[which.min(abs(rank(start_position) - median(rank(start_position))))]
  )
midpoints <- midpoints %>%
  left_join(k %>% dplyr::select(group, external_gene_name, start_position), 
            by = c("group", "mid_gene" = "external_gene_name")) %>%
  distinct(group, .keep_all = TRUE)

ggplot(scpar, aes(x = reorder(external_gene_name, start_position), y = ratio.xi.expression)) +
  geom_boxplot() +
  facet_grid(~group, scales = 'free_x', space = 'free') +
  theme_article() +
  geom_point(data = midpoints, 
             aes(x = reorder(mid_gene, start_position), y = mean_ae), 
             inherit.aes = FALSE, shape = 21, fill = "red", size = 3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

ggplot(scpar, aes(x = reorder(external_gene_name, start_position), y = ratio.xi.expression)) +
  geom_boxplot() +
  facet_grid(dataset ~ group, scales = 'free_x', space = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

ggplot(scpar, aes(x = reorder(external_gene_name, start_position), y = ratio.xi.expression)) +
  geom_boxplot() +
  facet_grid(dataset ~ group, scales = 'free_x', space = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

ggplot(scpar, aes(x = celltype.organ.sample, y = ratio.xi.expression)) +
  geom_boxplot() +
  facet_grid(.~group+external_gene_name, scales = 'free_x', space = 'free') +
  theme_article() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3),
        axis.title.x = element_blank())

mod = aov(data = scpar, ratio.xi.expression ~ group)
summary(mod)
TukeyHSD(mod)

# plot only averaged data

table(scpar$dataset)
pl = subset(scpar, dataset %in% c("scLinaX_AIDA", "scLinaX_multiome", "SanRoman_averaged","Tukanian_averaged"))

means <- pl %>%
  group_by(group) %>%
  summarise(mean_ae = mean(ratio.xi.expression, na.rm = TRUE))
midpoints <- pl %>%
  group_by(group) %>%
  summarise(
    mean_ae = mean(ratio.xi.expression, na.rm = TRUE),
    mid_gene = external_gene_name[which.min(abs(rank(start_position) - median(rank(start_position))))]
  )
midpoints <- midpoints %>%
  left_join(k %>% dplyr::select(group, external_gene_name, start_position), 
            by = c("group", "mid_gene" = "external_gene_name")) %>%
  distinct(group, .keep_all = TRUE)

ggplot(pl, aes(x = reorder(external_gene_name, start_position), y = ratio.xi.expression)) +
  geom_boxplot() +
  facet_grid(~group, scales = 'free_x', space = 'free') +
  theme_article() +
  geom_point(data = midpoints, 
             aes(x = reorder(mid_gene, start_position), y = mean_ae), 
             inherit.aes = FALSE, shape = 21, fill = "red", size = 3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

mod = aov(ratio.xi.expression ~ group, data = pl)
mod = aov(ratio.xi.expression ~ group + dataset + external_gene_name, data = pl)
summary(mod) 
TukeyHSD(mod) # A < B ; E < B C D padj<0.05

means <- pl %>%
  group_by(group, dataset) %>%
  summarise(mean_ae = mean(ratio.xi.expression, na.rm = TRUE))
midpoints <- pl %>%
  group_by(group, dataset) %>%
  summarise(
    mean_ae = mean(ratio.xi.expression, na.rm = TRUE),
    mid_gene = external_gene_name[which.min(abs(rank(start_position) - median(rank(start_position))))]
  )
midpoints <- midpoints %>%
  left_join(k %>% dplyr::select(group, external_gene_name, start_position), 
            by = c("group", "mid_gene" = "external_gene_name"))

ggplot(pl, aes(x = reorder(external_gene_name, start_position), y = ratio.xi.expression)) +
  geom_boxplot() +
  facet_grid(dataset~group, scales = 'free_x', space = 'free') +
  theme_article() +
  geom_point(data = midpoints, 
             aes(x = reorder(mid_gene, start_position), y = mean_ae), 
             inherit.aes = FALSE, shape = 21, fill = "red", size = 3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

# plot per sample data

table(scpar$dataset)
pl = subset(scpar, dataset %in% c("Tukiainen_persample", "SanRoman_persample"))

means <- pl %>%
  group_by(group) %>%
  summarise(mean_ae = mean(ratio.xi.expression, na.rm = TRUE))
midpoints <- pl %>%
  group_by(group) %>%
  summarise(
    mean_ae = mean(ratio.xi.expression, na.rm = TRUE),
    mid_gene = external_gene_name[which.min(abs(rank(start_position) - median(rank(start_position))))]
  )
midpoints <- midpoints %>%
  left_join(k %>% dplyr::select(group, external_gene_name, start_position), 
            by = c("group", "mid_gene" = "external_gene_name")) %>%
  distinct(group, .keep_all = TRUE)

ggplot(pl, aes(x = reorder(external_gene_name, start_position), y = ratio.xi.expression)) +
  geom_boxplot() +
  facet_grid(~group, scales = 'free_x', space = 'free') +
  theme_article() +
  geom_point(data = midpoints, 
             aes(x = reorder(mid_gene, start_position), y = mean_ae), 
             inherit.aes = FALSE, shape = 21, fill = "red", size = 3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

mod = aov(ratio.xi.expression ~ group, data = pl)
mod = aov(ratio.xi.expression ~ group + dataset + external_gene_name, data = pl)
summary(mod) 
TukeyHSD(mod) # A < B ; E < B C D padj<0.05

means <- pl %>%
  group_by(group, dataset) %>%
  summarise(mean_ae = mean(ratio.xi.expression, na.rm = TRUE))
midpoints <- pl %>%
  group_by(group, dataset) %>%
  summarise(
    mean_ae = mean(ratio.xi.expression, na.rm = TRUE),
    mid_gene = external_gene_name[which.min(abs(rank(start_position) - median(rank(start_position))))]
  )
midpoints <- midpoints %>%
  left_join(k %>% dplyr::select(group, external_gene_name, start_position), 
            by = c("group", "mid_gene" = "external_gene_name"))

ggplot(pl, aes(x = reorder(external_gene_name, start_position), y = ratio.xi.expression)) +
  geom_boxplot() +
  facet_grid(dataset~group, scales = 'free_x', space = 'free') +
  theme_article() +
  geom_point(data = midpoints, 
             aes(x = reorder(mid_gene, start_position), y = mean_ae), 
             inherit.aes = FALSE, shape = 21, fill = "red", size = 3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())


#### end ####


