##Function set names of ASVs


ASV_name_set <- function(ASVtab,TAXtab) 
  { # First put asv tab with sequence sin column headers
  out <- list()
  ASVnames=sprintf("ASV%s",seq(1:ncol(ASVtab))) ##### create a vector with new names
  sequences = colnames(ASVtab) #### extract the sequences
  sequencesnames=t(rbind(ASVnames, sequences)) #### create a data.frame with all this
  out$sequencesnames = sequencesnames
  colnames(ASVtab) = ASVnames 
  rownames(TAXtab) = ASVnames
  out$ASVtab = ASVtab
  out$TAXtab = TAXtab
  #### substitute the names in the ASV table
  return(out)
}




Plot_names  <- function(TAXtab) 
{ # First put asv tab with sequence sin column headers
  sequencesnames75 = TAXtab
  sequencesnames75$name = ifelse(!is.na(sequencesnames75$Species),paste(sequencesnames75$Genus,sequencesnames75$Species),sequencesnames75$Genus )
  sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$Family,sequencesnames75$name)
  sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$Order,sequencesnames75$name)
  sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$Class,sequencesnames75$name)
  sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$Phylum,sequencesnames75$name)
  sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$Kingdom,sequencesnames75$name)
  sequencesnames75$name = paste(rownames(sequencesnames75),sequencesnames75$name , sep = " ")
  #### substitute the names in the ASV table
  return(sequencesnames75$name)
}


generate.tax.summary.modified<-function(otu, tax.table) {
  tax.table$Species = ifelse(!is.na(tax.table$Species),paste(tax.table$Genus,tax.table$Species),paste(tax.table$Genus, "unknown"))
  tax.table[is.na(tax.table)] = "unknown"
  tax3.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Phylum, sum)))
  tax4.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Class, sum)))
  tax5.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Order, sum)))
  tax6.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Family, sum)))
  tax7.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Genus, sum)))
  tax8.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Species, sum)))
  return(list(tax3=tax3.summary, 
              tax4=tax4.summary, 
              tax5=tax5.summary, 
              tax6=tax6.summary, 
              tax7=tax7.summary, 
              tax8=tax8.summary))
}


#taxa.sum = generate.tax.summary.modified(rar.asv, as.data.frame(taxa))


`%nin%` = Negate(`%in%`)






kruskal_wallis_calc = function(ASVtab, vector) #ASVs as columns 
{
  out = data.frame()
  for (i in 1:ncol(ASVtab))
  {
    test = kruskal.test(ASVtab[,i] ~ vector)
    partial = c(colnames(ASVtab)[i], test$statistic, test$parameter,test$p.value)
    out = rbind(out, partial)
  }
  colnames(out) = c("gene","statistic","parameter","p_value")
  out$padj = p.adjust(out$p_value, method = "fdr")
  rownames(out) = out$gene
  out = out[,-1]
  return(out)
}




dunn_test_gabri = function(vector1, vector2)
  { 
              DT = FSA::dunnTest(vector1 ~ vector2,
              method="bh")      # Adjusts p-values for multiple comparisons;
              PT = DT$res
              return(rcompanion::cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05))
}

multiple_dunn_gabri = function(ASVtab, vector2)
{ 
                    df = data.frame()
                    df1 = data.frame()
                    for (i in 1:ncol(ASVtab))
                    {
                      DT = FSA::dunnTest(ASVtab[,i] ~ vector2,
                      method="bh")      # Adjusts p-values for multiple comparisons;
                      PT = DT$res
                      final = rcompanion::cldList(P.adj ~ Comparison,data = PT,threshold = 0.05)
                      final$Phylotype = colnames(ASVtab)[i]
                      PT$Phylotype = colnames(ASVtab)[i]
                      df = rbind(df, final)
                      df1 = rbind(df1,PT)
                    }
                    return(list(letters = df, tests = df1  ))
}




pearson_fuctions = function(ASVtab,vector)
{
  cor_results = data.frame()
  for ( i in 1:ncol(ASVtab))
  {
    test = cor.test( ASVtab[,i] ,  vector, method = "pearson")
    re_p = c(colnames(ASVtab)[i],test$estimate,  test$p.value)
    cor_results = rbind(cor_results,re_p)
  }
  colnames(cor_results) = c("ASV","estimate", "p.value")
  cor_results$q_value = p.adjust(cor_results$p.value, method = "fdr")
  return(cor_results)
}








spearman_fuctions = function(ASVtab,vector)
  {
    cor_results = data.frame()
    for ( i in 1:ncol(ASVtab))
    {
      test = cor.test( ASVtab[,i] ,  vector, method = "spearman")
      re_p = c(colnames(ASVtab)[i],test$estimate,  test$p.value)
      cor_results = rbind(cor_results,re_p)
    }
  colnames(cor_results) = c("ASV","estimate", "p.value")
  cor_results$q_value = p.adjust(cor_results$p.value, method = "fdr")
  return(cor_results)
}








wilcox_calc = function(ASVtab, vector) #ASVs as columns 
{
  out = data.frame()
  for (i in 1:ncol(ASVtab))
  {
    test = wilcox.test(ASVtab[,i] ~ vector, exact = FALSE)
    mean_flc = aggregate(ASVtab[,i], list(vector), FUN=mean, na.action = "na.omit") 
    l2cf = log2(mean_flc[1,2]/mean_flc[2,2])
    dem = mean_flc[1,1]
    partial = c(colnames(ASVtab)[i],  l2cf, dem,test$statistic, test$p.value)
    out = rbind(out, partial)
  }
  colnames(out) = c("gene", "l2cf", "numerator","statistic","p_value")
  out$padj = p.adjust(out$p_value, method = "fdr")
  out$l2cf = as.numeric(out$l2cf)
  rownames(out) = out$gene
  out = out[,-1]
  return(out)
}









anova_calc = function(ASVtab, vector) #ASVs as columns 
{
  out = data.frame()
  for (i in 1:ncol(ASVtab))
  {
    test = anova(lm(ASVtab[,i] ~ vector))
    partial = c(colnames(ASVtab)[i], test$Df[1], test$`F value`[1], test$`Pr(>F)`[1])
    out = rbind(out, partial)
  }
  colnames(out) = c("gene","df","F_value", "p_value")
  out$padj = p.adjust(out$p_value, method = "fdr")
  rownames(out) = out$gene
  out = out[,-1]
  return(out)
}





generate.tax.summary.modifiedGTDB<-function(otu, tax.table) {
  tax.table$Species = ifelse(is.na(tax.table$Species),paste(tax.table$Genus, "uknown", sep = " "),tax.table$Species)
  tax.table[is.na(tax.table)] = "unknown"
  tax3.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Phylum, sum)))
  tax4.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Class, sum)))
  tax5.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Order, sum)))
  tax6.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Family, sum)))
  tax7.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Genus, sum)))
  tax8.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$Species, sum)))
  return(list(tax3=tax3.summary, 
              tax4=tax4.summary, 
              tax5=tax5.summary, 
              tax6=tax6.summary, 
              tax7=tax7.summary, 
              tax8=tax8.summary))
}






Plot_names_GTDB  <- function(TAXtab) 
{ # First put asv tab with sequence sin column headers
  sequencesnames75 = TAXtab
  sequencesnames75$name = ifelse(!is.na(sequencesnames75$Species),sequencesnames75$Species,sequencesnames75$Genus )
  sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$Family,sequencesnames75$name)
  sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$Order,sequencesnames75$name)
  sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$Class,sequencesnames75$name)
  sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$Phylum,sequencesnames75$name)
  sequencesnames75$name = ifelse(is.na(sequencesnames75$name),sequencesnames75$Kingdom,sequencesnames75$name)
  sequencesnames75$name = paste(rownames(sequencesnames75),sequencesnames75$name , sep = " ")
  #### substitute the names in the ASV table
  return(sequencesnames75$name)
}










generate.tax.summary.sourmash<-function(otu, tax.table) {
  tax.table$Sspecies = ifelse(is.na(tax.table$species),paste(tax.table$Genus, "uknown", sep = " "),tax.table$species)
  tax.table[is.na(tax.table)] = "unknown"
  tax3.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$phylum, sum)))
  tax4.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$class, sum)))
  tax5.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$order, sum)))
  tax6.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$family, sum)))
  tax7.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$genus, sum)))
  tax8.summary<-as.data.frame(apply(otu, 1, function(x) by(x, tax.table$species, sum)))
  return(list(tax3=tax3.summary, 
              tax4=tax4.summary, 
              tax5=tax5.summary, 
              tax6=tax6.summary, 
              tax7=tax7.summary, 
              tax8=tax8.summary))
}





pearson_fuctions = function(ASVtab,vector)
{
  cor_results = data.frame()
  for ( i in 1:ncol(ASVtab))
  {
    test = cor.test( ASVtab[,i] ,  vector, method = "pearson")
    re_p = c(colnames(ASVtab)[i],test$estimate,  test$p.value)
    cor_results = rbind(cor_results,re_p)
  }
  colnames(cor_results) = c("ASV","estimate", "p.value")
  cor_results$q_value = p.adjust(cor_results$p.value, method = "fdr")
  return(cor_results)
}






anova_calc2 = function(ASVtab, vector) #ASVs as columns 
{
  out = data.frame()
  for (i in 1:ncol(ASVtab))
  {
    test = anova(lm(vector ~ ASVtab[,i]))
    partial = c(colnames(ASVtab)[i], test$Df[1], test$`F value`[1], test$`Pr(>F)`[1])
    out = rbind(out, partial)
  }
  colnames(out) = c("gene","df","F_value", "p_value")
  out$padj = p.adjust(out$p_value, method = "fdr")
  rownames(out) = out$gene
  out = out[,-1]
  return(out)
}




## Multiple permanovas -----
permanova_multi = function(ASVtab, metadata) #ASVs as columns 
{
  out = data.frame()
  for (i in 1:ncol(metadata))
  {
    a = adonis2(ASVtab ~ metadata[,i], na.action = na.omit)
    partial = c(colnames(metadata)[i],a$Df[1], a$Df[2], a$F[1], a$R2[1], a$`Pr(>F)`[1])
    out = rbind(out, partial)
  }
  colnames(out) = c("Variable","Df","DF_res", "F.Model", "R2", "p_value")
  out$padj = p.adjust(out$p_value, method = "fdr")
  rownames(out) = out$gene
  return(out)
}




select_top = function(ASVtab, number, mean_or_sum)
{
  if(mean_or_sum == "mean")
  {
    sums = colMeans(ASVtab)
    top = sort(sums, decreasing = TRUE)[1:number]
  }
  else
  {
    sums = colSums(ASVtab)
    top = sort(sums, decreasing = TRUE)[1:number]
  }
  ASVtab = ASVtab[,colnames(ASVtab) %in% names(top)]
  return(ASVtab)
}




Css_trn = function(ASVtab)
{  
  if (!require(metagenomeSeq)) {stop("metagenomeSeq not installed")}
  metaSeqObject    = newMRexperiment(t(ASVtab) )
  metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
  ASVtab_16s_CSS = data.frame(t(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE)))
  return(ASVtab_16s_CSS)
}


perform_permanova <- function(data, asvtab) {
  results_df <- data.frame(Variable = character(0), P_value = numeric(0))
  distance <- asvtab
  for (l in 1:ncol(data)) {
                variable <- data[, l]
                data_c <- cbind(distance, variable)
                colnames(data_c)[ncol(data_c)] <- colnames(data)[l]
                data_c <- na.omit(data_c)
  # Splitting the data_c into distance matrix and variable
distance1 <- data_c[, 1:ncol(distance)]
# Create the formula
formula <- as.formula(paste("distance1 ~", colnames(data_c)[ncol(data_c)]))
                             
# Perform PERMANOVA
permanova_result <- adonis(formula, data = data_c)
 p_value <- permanova_result$aov.tab$"Pr(>F)"[1]
results_df <- rbind(results_df, data.frame(Variable = colnames(data)[l], P_value = p_value))
}
                           
return(results_df)
}















Order_ASV_mean_sum = function(vector_asv, ASVtab, mean_or_sum)
{
  ASVtab = ASVtab[,colnames(ASVtab) %in% vector_asv]
  if(mean_or_sum == "mean")
  {
    sums = colMeans(ASVtab)
    top = sort(sums, decreasing = TRUE)
  }
  else
  {
    sums = colSums(ASVtab)
    top = sort(sums, decreasing = TRUE)
  }
  ASVtab = ASVtab[,colnames(ASVtab) %in% names(top)]
  ASVtab = ASVtab[,match(names(top),colnames(ASVtab))]
  return(ASVtab)
}  






       