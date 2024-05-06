library("tidyverse")
library("ggpubr")
library("rstatix")

args<-c("Path_to_replicates", "pattern")

list_files<-list.files(args[1], pattern = as.character(args[2]), full.names = TRUE)

df<-data.frame()
Model<-c()
AICc<-c()
Replicate<-c()
for( i in list_files){
        iteration<-gsub(paste(args[1], "/df_", as.character(args), sep =""),'', i)
        iteration<-gsub('.tsv', '',  iteration)
        file<-read.csv(i, sep ="\t")
        for (j in 1:nrow(file)){
            Model<-c(Model, rownames(file[j,]))
            AICc<-c(AICc, file[j,2])   
            Replicate<-c(Replicate, as.numeric(iteration))
        }
}
df<-as.data.frame(cbind(Model, AICc, Replicate))
df$AICc<-as.numeric(df$AICc)

tab_test_shap <- df %>%
  group_by(Model) %>%
  shapiro_test(AICc)

tab_test_shap

tab_test_mean <-df %>%
  group_by(Model) %>%
  get_summary_stats(AICc, type = "mean_sd")
tab_test_mean

df_dups <- df[c("Model", "Replicate")]

df<-df[!duplicated(df_dups),]

res.aov <- df %>% friedman_test(AICc ~ Model |Replicate)
    pwc <- df %>%
  wilcox_test(AICc ~ Model, paired = TRUE, p.adjust.method = "bonferroni")

res.aov
pwc

pwc2 <- pwc %>% add_xy_position(x = "time")
    ggboxplot(df, x = "Model", y = "AICc", add = "point") +
      stat_pvalue_manual(pwc2, hide.ns = TRUE) +
      labs(
        subtitle = get_test_label(res.aov,  detailed = TRUE),
        caption = get_pwc_label(pwc2)
          )
ggsave("result_AICc_boxplot.pdf")
