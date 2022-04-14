library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
remotes::install_github("teunbrand/ggh4x")

###########################################################################
########################## plot TE master library #########################
###########################################################################
annotation <- read.gff("~/Github/jingxuan/data/te_lib/full_length_from_genome.sorted.gff", na.strings = c(".", "?"), GFF3 = TRUE)
# get_plot <- function(x){
#   setcolors <- c(
#     "TY1"="#A6CEE3", "Ty1_mosaic"="#1F78B4", "Ty1_prime"="#B2DF8A", "Ty1p_nw"="#33A02C", "Ty1p_ow"="#FB9A99", "Ty2"="#E31A1C",
#     "Ty3"="#FDBF6F", "Ty3p_ow"="#FF7F00", "Ty4"="#CAB2D6", "Ty4p_ow"="#6A3D9A", "Tsu4p_nw"="#E1BE6A", "Ty5"="#B15928", "Ty5p_ow"="#000000"
#   )
#   df <- annotation %>% filter(seqid==x) %>% mutate(name = case_when(attributes=="Name=GAG" ~ "GAG", attributes=="Name=POL" ~ "POL", attributes=="Name=GAG-POL" ~ "GAG-POL", attributes=="Name=5' long terminal repeat" ~ "LTR", attributes=="Name=3' long terminal repeat" ~ "LTR", TRUE ~ attributes))
#   df$name <- factor(df$name, levels = c("POL","GAG","GAG-POL","LTR"))
#   length <- filter(df, type=="LTR_retrotransposon")$end %>% as.numeric()
#   p <- ggplot(df) +
#     geom_segment(data=filter(df, !type=="LTR_retrotransposon"), aes(x=start, xend=end, y=name, yend=name, color = seqid), size=3) +
#     theme_minimal() +
#     theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = "none", axis.line.x = element_line(), axis.ticks.x = element_line()) +
#     xlab(NULL) + ylab(NULL) +
#     scale_x_continuous(breaks=c(1,3000, length)) +
#     scale_color_manual(values = setcolors) +
#     ggtitle(x)
#   return(p)
# }
#   
# p_TY1 <- get_plot("TY1")
# p_TY1canonical <- get_plot("Ty1_mosaic")
# p_TY1prime <- get_plot("Ty1_prime")
# p_TY1p_nw <- get_plot("Ty1p_nw")
# p_TY1p_ow <- get_plot("Ty1p_ow")
# p_TY2 <- get_plot("Ty2")
# p_TY3 <- get_plot("Ty3")
# p_TY3_1p <- get_plot("Ty3p_ow")
# p_TY4 <- get_plot("Ty4")
# p_TY4p_ow <- get_plot("Ty4p_ow")
# p_TSU4 <- get_plot("Tsu4p_nw")
# p_TY5c <- get_plot("Ty5")
# p_TY5p <- get_plot("Ty5p_ow")
# 
# ggpubr::ggarrange(p_TY1, p_TY1canonical, p_TY1prime, p_TY1p_nw, p_TY1p_ow, p_TY2, p_TY3, p_TY3_1p, p_TY4, p_TY4p_ow, p_TSU4, p_TY5c, p_TY5p, ncol = 1, nrow = 13)
# ggsave("~/Github/jingxuan/data/te_lib/master_lib.jpg", device = "jpg", width = 15, height = 45, units = "cm", dpi = 320)
# ggsave("~/Github/jingxuan/data/te_lib/master_lib.pdf", device = "pdf", width = 15, height = 45, units = "cm", dpi = 320)


df <- annotation %>% 
  mutate(name = case_when(attributes=="Name=GAG" ~ "GAG", attributes=="Name=POL" ~ "POL", attributes=="Name=GAG-POL" ~ "GAG-POL", attributes=="Name=5' long terminal repeat" ~ "LTR", attributes=="Name=3' long terminal repeat" ~ "LTR", TRUE ~ attributes)) %>% 
  separate(attributes, c("family","organism","strain","coordinates"), remove = F, sep = ";", fill = "right") %>%
  mutate(family = sub("^.*=","",family,perl = T), organism = sub("^.*=","",organism,perl = T), strain = sub("^.*=","",strain,perl = T), coordinates = sub("^.*=","",coordinates,perl = T)) %>%
  mutate(title = paste0(family,", ", sub("Saccharomyces","S.",organism), ", ", strain, "\n", coordinates))
# order of annotations
df$name <- factor(df$name, levels = c("POL","GAG","GAG-POL","LTR"))
# order of ty family
df$seqid <- factor(df$seqid, levels = c("Ty1_canonical","Ty1_mosaic","Ty1_prime","Ty1p_nw","Ty1p_ow","Ty2","Ty3","Ty3p_ow","Ty4","Ty4p_ow","Tsu4p_nw","Ty5","Ty5p_ow"))
# colors for ty family
setcolors <- c(
  "Ty1_canonical"="#A6CEE3", "Ty1_mosaic"="#1F78B4", "Ty1_prime"="#B2DF8A", "Ty1p_nw"="#33A02C", "Ty1p_ow"="#FB9A99", "Ty2"="#E31A1C",
  "Ty3"="#FDBF6F", "Ty3p_ow"="#FF7F00", "Ty4"="#CAB2D6", "Ty4p_ow"="#6A3D9A", "Tsu4p_nw"="#E1BE6A", "Ty5"="#B15928", "Ty5p_ow"="#000000"
)
# modify title to include more details
family <- filter(df, type=="LTR_retrotransposon")$seqid; title <- filter(df, type=="LTR_retrotransposon")$title
names(title) <- family
panel_labels <- labeller(seqid = title)
# modify scale for each ty family
xscales <- list(
  seqid == "Ty1_canonical" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty1_canonical")$end))),
  seqid == "Ty1_mosaic" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty1_mosaic")$end))),
  seqid == "Ty1_prime" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty1_prime")$end))),
  seqid == "Ty1p_nw" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty1p_nw")$end))),
  seqid == "Ty1p_ow" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty1p_ow")$end))),
  seqid == "Ty2" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty2")$end))),
  seqid == "Ty3" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty3")$end))),
  seqid == "Ty3p_ow" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty3p_ow")$end))),
  seqid == "Ty4" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty4")$end))),
  seqid == "Ty4p_ow" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty4p_ow")$end))),
  seqid == "Tsu4p_nw" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Tsu4p_nw")$end))),
  seqid == "Ty5" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty5")$end))),
  seqid == "Ty5p_ow" ~ scale_x_continuous(breaks = c(1,3000,as.numeric(filter(df, type=="LTR_retrotransposon"&seqid=="Ty5p_ow")$end)))
)
# generate plot
ggplot(df) +
  geom_segment(data=filter(df, !type=="LTR_retrotransposon"), aes(x=start, xend=end, y=name, yend=name, color = seqid), size=3) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        legend.position = "none", 
        axis.line.x = element_line(), axis.ticks.x = element_line(),
        strip.text.x = element_text(face = "bold")) +
  xlab(NULL) + ylab(NULL) +
  # scale_x_continuous(breaks=c(1,3000, length)) +
  scale_color_manual(values = setcolors) +
  facet_wrap(~ seqid, scales = "free", ncol = 1, labeller = panel_labels) +
  ggh4x::facetted_pos_scales(x = xscales)
# save files
ggsave("~/Github/jingxuan/data/te_lib/master_lib.jpg", device = "jpg", width = 15, height = 45, units = "cm", dpi = 320)
ggsave("~/Github/jingxuan/data/te_lib/master_lib.pdf", device = "pdf", width = 15, height = 45, units = "cm", dpi = 320)

