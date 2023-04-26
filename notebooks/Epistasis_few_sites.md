Epistasis analysis
================
Santiago Herrera
2023-04-25

Functions
---------

``` r
# imports a fasta-formatted alignment and converts it into a matrix
import_alignment <- function(dir) {
  
  fasta <- readLines(dir)
  
  taxa <- fasta[2L * (1L:(length(fasta) / 2L)) - 1L]
  taxa <- sub('>', '', taxa)
  
  sequences <- fasta[2L * (1L:(length(fasta) / 2L))]
  sequences <- t(vapply(strsplit(sequences, split = ''), identity, character(nchar(sequences[1L]))))
  rownames(sequences) <- taxa
  
  sequences
}

# Filter an alignment matrix to include only certain taxa and sites
filter_alignment_by_taxa_and_sites <- function(taxa,sites,aln){
  m <- aln[row.names(aln) %in% taxa,]
  s <- sapply(sites, function(x) paste0("V",x))
  m2 <- m[,colnames(m) %in% s]
  return(m2)
}

# For two taxa get the number (and percentage) of amino acid differences
get_pairwise_dist <- function(aln,df){
  aa_dist <- c()
  perc_dist <- c()
  for(i in 1:dim(df)[1]){
    t <- c(df[i,3],df[i,4])
    a <- aln[row.names(aln) %in% t,]
    diffs <- 0
    for(j in 1:dim(a)[2]){
      if(a[1,j] != a[2,j]) diffs <- diffs + 1
    }
    aa_dist <- c(aa_dist,diffs)
    perc_dist <- c(perc_dist,((diffs/dim(a)[2])*100))
  }
  return(list(aa_dist,perc_dist))
}

# Creates a matrix with recoded states as 0 or 1 depending on whether a pairwise construct comparison have the same or different amino acid state at that site, respectively.
linear_model_matrix <- function(aln,df){
  
  matrix_lm <- matrix(data = NA,nrow = dim(df)[1],ncol = dim(aln)[2]+2)
  effects <- df[,8] # Y values for ANOVA
  covariate <- df[,10] # covariate to correct ANOVA
  matrix_lm[,1] <- effects
  matrix_lm[,2] <- covariate
  
  # For every site in alignment:
  for(i in 1:dim(aln)[2]){
    # If site has no variation in alignment, continue to next site:
    if(length(unique(aln[,i])) == 1) next
    
    # If more than 50% of taxa have gaps at that site, continue to next site:
    if(sum(aln[,1]=="-")/length(aln[,i]) > 0.5) next
    
    aa_states <- c()
    
    # For every pairwise construct comparison:
    for(j in 1:dim(df)[1]){
      t <- c(df[j,3],df[j,4])
      a <- aln[row.names(aln) %in% t, i]  # sub-alignment: compare single variable site amongst two taxa
      # Is the amino acid state the same? --> Recode states as binary (0 or same, and 1 or diff.) to increase statistical power
      if(a[1] == a[2])  aa_states <- c(aa_states,0)
      if(a[1] != a[2])  aa_states <- c(aa_states,1)
    }
    
    # Check points to ensure there's enough power and observations to perform an ANOVA:
      # If site has the same recoded state amongst focal taxa, continue (either all states are the same, or all states are different).
      # If only one observation in either recoded binary state, continue.
    if(length(unique(aa_states)) == 1) next
    if(table(aa_states)[1] < 2) next
    if(table(aa_states)[2] < 2) next
    matrix_lm[,i+2] <- as.numeric(aa_states) # save recoded states to matrix
  }
  matrix_lm <- as.data.frame(matrix_lm)
  colnames(matrix_lm) <- c("effect","covariate",sapply(seq(1,dim(aln)[2],1), function(x) paste0("V",x)))
  
  matrix_lm <- matrix_lm[ , colSums(is.na(matrix_lm)) == 0]
  return(matrix_lm)
}

# This function performs an ANOVA for every variant site between the functional constructs.
anova_per_site <- function(mat){
  # mat: a matrix with amino states recoded. (output of `joint_linear_model_matrix` function)
  
  position <- c()
  effects <- mat[,1] # Y values for ANOVA
  covariate <- mat[,2] # covariate to correct ANOVA
  pvals <- c()
  Fvals <- c()
  R2 <- c()
  
  for(i in 3:dim(mat)[2]){
    # Do ANOVA at site i
    aa_states <- as.factor(mat[,i])
    m <- aov(effects~covariate+aa_states)
    p <- summary(m)[[1]][["Pr(>F)"]][2]
    f <- summary(m)[[1]][["F value"]][2]
    r2 <- summary(m)[[1]][["Sum Sq"]][2]/(summary(m)[[1]][["Sum Sq"]][1]+summary(m)[[1]][["Sum Sq"]][2]+summary(m)[[1]][["Sum Sq"]][3])
    position <- c(position,colnames(mat)[i])
    pvals <- c(pvals,p)
    Fvals <- c(Fvals,f)
    R2 <- c(R2,r2)
  }
  df <- data.frame(position,Fvals,R2,pvals)
  df$site <- match_site(df)
  return(df)
}

#Function to match alignment position to sheep A1 site number (see "A1_AA_alignment_convergence_modified_PA_codeml.fasta")
##  Blocks of alignment PA trimmed      Complete alignment   
##  1                   1 - 15          14 - 28               1 : 1
##  2                   16 - 24         29 - 47               Not 1 : 1
##  3                   25 - 114        48 - 137              1 : 1
##  4                   115 - end       139 - end             1 : 1

match_site <- function(table){
  sites <- c()
  for(i in 1:dim(table)[1]){
    position <- as.numeric(table[i,1] %>% gsub("V","",.))
    site <- NA
    if(position >= 14 && position <= 28) site <- position - 13
    if(position >= 48 && position <= 137) site <- position - 23
    if(position >= 139) site <- position - 24
    sites <- c(sites,site)
  }
  return(sites)
}

# Extract R^2 from a linear model
get_R2 <- function(model){
  ss_cov <- summary(model)[[1]][["Sum Sq"]][1]
  ss_res <- tail(summary(model)[[1]][["Sum Sq"]],1)
  ss_factors <- summary(model)[[1]][["Sum Sq"]][-c(which(summary(model)[[1]][["Sum Sq"]]==ss_cov),which(summary(model)[[1]][["Sum Sq"]]==ss_res))]
  R2 <- sum(ss_factors)/sum(summary(model)[[1]][["Sum Sq"]]) 
  R2
}

# Extract adjusted R^2 from a linear model (Adjusted R2 also indicates how well terms fit a curve or line, but adjusts for the number of terms in a model.)
adj_R2 <- function(model){
  n <- length(model$model$effect)
  k <- model$rank-1
  adjr2 <- 1-(((1-get_R2(model))*(n-1))/(n-k-1))
  adjr2
} 
```

This script has the analyses performed in fig 5. pertaining the finding that a small number of sites account for a large proportion of the differences in pleiotropic effects of the same substitution on divergent ATP1A1 backgrounds.

Import data
-----------

``` r
# Import alignment of extant sequences
alignment <- "../data/alignments/A1-A3_AA_alignment_convergence_modified_codeml.fasta"
aln <- import_alignment(alignment)
colnames(aln) <- sapply(seq(1,1040,1), function(x) paste0("V",x))

# Import dataset with the pairwise construct comparisons (Distance and Perc_dist refer to the full protein)
df_effect <- read.csv("../data/mut_effects/Percent_change_shared_derived_state.csv",header = T)

# Names of taxa with functional information (use to make sub-alignmentes)
taxa <- c("bird|Struthio_camelus_A1","frog|Leptodactylus_macrosternumS_A1","mammal|Rattus_norvegicus_A1",
          "snake|Xenodon_rhabdocephalus_A1","lizard|Tupinambis_teguxin_A1","lizard|Varanus_exanthemathicus_A1",
          "snake|Rhabdophis_subminiatus_A1","mammal|Chinchilla_lanigera_A1")
```

The effect of divergence at all sites
-------------------------------------

Relationship between percent change (when change to the same derived amino acid state) and number of amino acid differences for the entire protein.

``` r
# Statistical analyses:
# Correlation between percent change and number of amino acid differences for the entire protein 
cor.test(df_effect$Distance,df_effect$PercentDiff)
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  df_effect$Distance and df_effect$PercentDiff
    ## t = -1.3052, df = 9, p-value = 0.2242
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.8059457  0.2641394
    ## sample estimates:
    ##        cor 
    ## -0.3989494

``` r
# Plot
my_pallette <- brewer.pal(6,"Spectral")

ggplot(df_effect,aes(x=Distance,y=PercentDiff),color=AAState) + 
  geom_point(size=4,shape=21,aes(fill = AAState)) +
  scale_fill_manual(values = my_pallette) +
  labs(x="# Amino acid differences between ATP1A1 orthologs",y="Difference in effect on two backgrounds (%)",fill="AA State") +
  theme(
    panel.background = element_rect(fill = "white",colour = "white"),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid.major = element_line(),
    axis.text.x = element_text(face="bold", size=10),
    axis.text.y = element_text(face="bold", size=10), 
    axis.title = element_text(face="bold", size=13),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size=16),
    legend.title = element_text(size=16)) +
  annotate("text",x=65,y=180,label=expression(paste(rho, " = -0.39, P > 0.05")))
```

![](Epistasis_few_sites_files/figure-markdown_github/allsiteseffect-1.png)

ANOVA per site
--------------

Here we will perform an ANOVA for every variable site in the wild-type sequence amongst the constructs with functional data.

``` r
# Add covariate to functional data frame: which site is being mutated
df_effect$covariate <- gsub("[A-Z]","",df_effect$AAState)

# Make sub-alignment (Taxa with functional data and full ATP1A1 sequence)
constructs_aln <- filter_alignment_by_taxa_and_sites(taxa,seq(1,dim(aln)[2],1),aln)

# Build recoded matrix
matrix_anova <- linear_model_matrix(constructs_aln,df_effect) 

# Remove polymorphic positions within the first 68 sites of the alignment (poor alignment quality) and sites 111 and 122.
sites_to_remove <- c("V1","V2","V9","V10","V12","V14","V15","V16","V17","V19","V20","V21","V22","V23","V25","V26","V27","V28","V29","V30","V31","V32","V33","V34","V35","V36","V43","V44","V45","V46","V48","V56","V57","V62","V63","V69","V76","V80","V134","V146")
matrix_anova <- matrix_anova[, !(colnames(matrix_anova) %in% c(sites_to_remove))]

# Perform and ANOVA per site 
df_anova <- anova_per_site(matrix_anova) 
df_anova$LogP <- -log(df_anova$pvals)
df_anova <- df_anova[!(is.na(df_anova$site)),] # remove sites that have unclear match to reference sequence
```

Correlation amongst sites
-------------------------

Because we have relatively few construct comparisons relative to the number of variable sites, some sites have strong correlations between constructs (i.e., they show the same pattern of variation, therefore have the same signal).

``` r
# Build correlation matrix amogst sites based on their variation pattern amongst constructs
matrix_sites_only <- matrix_anova[,!(colnames(matrix_anova) %in% c("effect","covariate"))]
matrix_sites_only[] <- sapply(matrix_sites_only,as.numeric)
cor_matrix <- cor(matrix_sites_only)
cor_matrix[!lower.tri(cor_matrix)] <- 0 # Set upper diagonal to zero

# Plot absolute correlation between sites to show clusters of highly colinear sites
pheatmap(abs(cor(matrix_sites_only)),cluster_rows = T, border_color="black",fontsize_row=3,fontsize_col=3,cellheight=2,main = "Clusters of polymorphic sites")
```

![](Epistasis_few_sites_files/figure-markdown_github/corrsites-1.png)

``` r
# Remove highly colinear variables (sites with an absolute correlation > 0.8)
matrix_anova_indep_1 <- matrix_sites_only[, !apply(cor_matrix, 2, function(x) any(abs(x) > 0.8, na.rm = TRUE))] # Selects an arbritary site as representatives
matrix_anova_indep_1$effect <- matrix_anova$effect
matrix_anova_indep_1$covariate <- matrix_anova$covariate

# Remove highly colinear variables (sites with an absolute correlation > 0.99)
matrix_anova_indep_2 <- matrix_sites_only[, !apply(cor_matrix, 2, function(x) any(abs(x) > 0.99, na.rm = TRUE))]
matrix_anova_indep_2$effect <- matrix_anova$effect
matrix_anova_indep_2$covariate <- matrix_anova$covariate

# Sort ANOVA results by the marginal R^2
df_anova <- arrange(df_anova,desc(R2))

# Add group category to sites: grouped by the cluster they belong to
group_sites1 <- matrix_anova_indep_1 %>% select(-c(effect,covariate)) %>% with(colnames(.))
group_sites1 <- group_sites1[order(match(group_sites1,df_anova$position[which(df_anova$position %in% group_sites1)]))] # order representative sites based on R2
group1 <- seq(1,length(group_sites1),1)
group_sites1 <- data.frame(position=group_sites1,group80=group1)

group_sites2 <- matrix_anova_indep_2 %>% select(-c(effect,covariate)) %>% with(colnames(.))
group_sites2 <- group_sites2[order(match(group_sites2,df_anova$position[which(df_anova$position %in% group_sites2)]))] # order representative sites based on R2
group2 <- seq(1,length(group_sites2),1)
group_sites2 <- data.frame(position=group_sites2,group99=group2)

df_anova <- df_anova %>% full_join(.,group_sites1,by="position") %>% full_join(.,group_sites2,by="position") %>% fill(group80,.direction="up") %>% fill(group99,.direction="up") %>% mutate(group80 = replace_na(group80,24))

# Show ANOVA table with assigned groups
knitr::kable(df_anova)
```

| position |       Fvals|         R2|      pvals|  site|       LogP|  group80|  group99|
|:---------|-----------:|----------:|----------:|-----:|----------:|--------:|--------:|
| V250     |  15.7110034|  0.6626033|  0.0041554|   226|  5.4833534|        1|        1|
| V295     |  15.7110034|  0.6626033|  0.0041554|   271|  5.4833534|        1|        1|
| V314     |  15.7110034|  0.6626033|  0.0041554|   290|  5.4833534|        1|        1|
| V430     |  15.7110034|  0.6626033|  0.0041554|   406|  5.4833534|        1|        1|
| V669     |  15.7110034|  0.6626033|  0.0041554|   645|  5.4833534|        1|        1|
| V856     |  15.7110034|  0.6626033|  0.0041554|   832|  5.4833534|        1|        1|
| V125     |  11.9761117|  0.5995211|  0.0085596|   102|  4.7607045|        1|        2|
| V203     |  11.9761117|  0.5995211|  0.0085596|   179|  4.7607045|        1|        2|
| V248     |  11.9761117|  0.5995211|  0.0085596|   224|  4.7607045|        1|        2|
| V251     |  11.9761117|  0.5995211|  0.0085596|   227|  4.7607045|        1|        2|
| V440     |  11.9761117|  0.5995211|  0.0085596|   416|  4.7607045|        1|        2|
| V545     |  11.9761117|  0.5995211|  0.0085596|   521|  4.7607045|        1|        2|
| V912     |  11.9761117|  0.5995211|  0.0085596|   888|  4.7607045|        1|        2|
| V1003    |  11.9761117|  0.5995211|  0.0085596|   979|  4.7607045|        1|        2|
| V585     |   8.8945423|  0.5264738|  0.0175368|   561|  4.0434523|        2|        3|
| V303     |   7.7496571|  0.4920520|  0.0237839|   279|  3.7387446|        2|        4|
| V515     |   3.9499791|  0.3305425|  0.0820958|   491|  2.4998684|        3|        5|
| V143     |   3.6622381|  0.3140250|  0.0920053|   119|  2.3859089|        3|        6|
| V201     |   3.6622381|  0.3140250|  0.0920053|   177|  2.3859089|        3|        6|
| V457     |   3.6622381|  0.3140250|  0.0920053|   433|  2.3859089|        3|        6|
| V599     |   3.6622381|  0.3140250|  0.0920053|   575|  2.3859089|        3|        6|
| V899     |   3.6622381|  0.3140250|  0.0920053|   875|  2.3859089|        3|        6|
| V492     |   3.4820092|  0.3032576|  0.0990116|   468|  2.3125182|        4|        7|
| V439     |   3.0636900|  0.2769137|  0.1181716|   415|  2.1356171|        5|        8|
| V277     |   2.2878510|  0.2223836|  0.1688439|   253|  1.7787804|        5|        9|
| V433     |   2.2878510|  0.2223836|  0.1688439|   409|  1.7787804|        5|        9|
| V485     |   2.2878510|  0.2223836|  0.1688439|   461|  1.7787804|        5|        9|
| V682     |   2.2878510|  0.2223836|  0.1688439|   658|  1.7787804|        5|       10|
| V1029    |   2.2878510|  0.2223836|  0.1688439|  1005|  1.7787804|        5|       10|
| V436     |   1.9218622|  0.1936996|  0.2030625|   412|  1.5942414|        6|       11|
| V516     |   1.8740797|  0.1897977|  0.2082085|   492|  1.5692151|        6|       12|
| V139     |   1.6181335|  0.1682376|  0.2390806|   115|  1.4309544|        6|       13|
| V331     |   1.4285994|  0.1515175|  0.2662197|   307|  1.3234332|        7|       14|
| V276     |   1.3659169|  0.1458390|  0.2761590|   252|  1.2867786|        7|       15|
| V426     |   1.3659169|  0.1458390|  0.2761590|   402|  1.2867786|        7|       15|
| V491     |   1.3659169|  0.1458390|  0.2761590|   467|  1.2867786|        7|       15|
| V494     |   1.3659169|  0.1458390|  0.2761590|   470|  1.2867786|        7|       15|
| V518     |   1.3659169|  0.1458390|  0.2761590|   494|  1.2867786|        7|       15|
| V520     |   1.3659169|  0.1458390|  0.2761590|   496|  1.2867786|        7|       15|
| V898     |   1.3659169|  0.1458390|  0.2761590|   874|  1.2867786|        7|       15|
| V938     |   1.3271582|  0.1422895|  0.2825731|   914|  1.2638179|        7|       16|
| V700     |   1.3271582|  0.1422895|  0.2825731|   676|  1.2638179|        7|       17|
| V540     |   1.3201476|  0.1416444|  0.2837562|   516|  1.2596399|        7|       18|
| V580     |   1.3039855|  0.1401533|  0.2865110|   556|  1.2499783|        8|       19|
| V597     |   1.0485167|  0.1158771|  0.3358121|   573|  1.0912036|        8|       20|
| V298     |   0.9852040|  0.1096473|  0.3499863|   274|  1.0498613|        9|       21|
| V490     |   0.9677981|  0.1079192|  0.3540429|   466|  1.0383371|        9|       22|
| V204     |   0.9677981|  0.1079192|  0.3540429|   180|  1.0383371|        9|       23|
| V453     |   0.9677981|  0.1079192|  0.3540429|   429|  1.0383371|        9|       23|
| V463     |   0.9677981|  0.1079192|  0.3540429|   439|  1.0383371|        9|       23|
| V495     |   0.9677981|  0.1079192|  0.3540429|   471|  1.0383371|        9|       23|
| V509     |   0.9677981|  0.1079192|  0.3540429|   485|  1.0383371|        9|       23|
| V541     |   0.9677981|  0.1079192|  0.3540429|   517|  1.0383371|        9|       23|
| V578     |   0.9677981|  0.1079192|  0.3540429|   554|  1.0383371|        9|       23|
| V591     |   0.9677981|  0.1079192|  0.3540429|   567|  1.0383371|        9|       23|
| V693     |   0.9677981|  0.1079192|  0.3540429|   669|  1.0383371|        9|       23|
| V695     |   0.9677981|  0.1079192|  0.3540429|   671|  1.0383371|        9|       23|
| V749     |   0.9677981|  0.1079192|  0.3540429|   725|  1.0383371|        9|       23|
| V945     |   0.9677981|  0.1079192|  0.3540429|   921|  1.0383371|        9|       23|
| V1017    |   0.9677981|  0.1079192|  0.3540429|   993|  1.0383371|        9|       23|
| V311     |   0.9247202|  0.1036133|  0.3643990|   287|  1.0095057|        9|       24|
| V131     |   0.8089222|  0.0918298|  0.3947044|   108|  0.9296181|       10|       25|
| V189     |   0.8089222|  0.0918298|  0.3947044|   165|  0.9296181|       10|       25|
| V1036    |   0.8089222|  0.0918298|  0.3947044|  1012|  0.9296181|       10|       25|
| V141     |   0.7653443|  0.0873147|  0.4071577|   117|  0.8985548|       11|       26|
| V598     |   0.6830762|  0.0786675|  0.4324863|   574|  0.8382046|       11|       27|
| V124     |   0.6725032|  0.0775442|  0.4359306|   101|  0.8302723|       11|       28|
| V142     |   0.6725032|  0.0775442|  0.4359306|   118|  0.8302723|       11|       28|
| V289     |   0.6725032|  0.0775442|  0.4359306|   265|  0.8302723|       11|       28|
| V118     |   0.6016437|  0.0699451|  0.4602641|    95|  0.7759549|       12|       29|
| V551     |   0.6016437|  0.0699451|  0.4602641|   527|  0.7759549|       12|       29|
| V95      |   0.5721575|  0.0667460|  0.4710885|    72|  0.7527094|       12|       30|
| V136     |   0.5721575|  0.0667460|  0.4710885|   113|  0.7527094|       12|       30|
| V176     |   0.5721575|  0.0667460|  0.4710885|   152|  0.7527094|       12|       30|
| V489     |   0.5721575|  0.0667460|  0.4710885|   465|  0.7527094|       12|       30|
| V536     |   0.5721575|  0.0667460|  0.4710885|   512|  0.7527094|       12|       30|
| V692     |   0.5721575|  0.0667460|  0.4710885|   668|  0.7527094|       12|       30|
| V903     |   0.5721575|  0.0667460|  0.4710885|   879|  0.7527094|       12|       30|
| V915     |   0.5721575|  0.0667460|  0.4710885|   891|  0.7527094|       12|       30|
| V431     |   0.5721575|  0.0667460|  0.4710885|   407|  0.7527094|       13|       31|
| V1021    |   0.5721575|  0.0667460|  0.4710885|   997|  0.7527094|       13|       31|
| V226     |   0.5218396|  0.0612355|  0.4906328|   202|  0.7120594|       14|       32|
| V513     |   0.5218396|  0.0612355|  0.4906328|   489|  0.7120594|       14|       32|
| V537     |   0.5218396|  0.0612355|  0.4906328|   513|  0.7120594|       14|       32|
| V992     |   0.5218396|  0.0612355|  0.4906328|   968|  0.7120594|       14|       32|
| V486     |   0.4326106|  0.0513021|  0.5291777|   462|  0.6364310|       15|       33|
| V691     |   0.3980993|  0.0474035|  0.5456719|   667|  0.6057373|       16|       34|
| V479     |   0.3547464|  0.0424604|  0.5679017|   455|  0.5658069|       17|       35|
| V519     |   0.3547464|  0.0424604|  0.5679017|   495|  0.5658069|       17|       35|
| V694     |   0.3547464|  0.0424604|  0.5679017|   670|  0.5658069|       17|       35|
| V701     |   0.3547464|  0.0424604|  0.5679017|   677|  0.5658069|       17|       35|
| V754     |   0.3547464|  0.0424604|  0.5679017|   730|  0.5658069|       17|       35|
| V905     |   0.3547464|  0.0424604|  0.5679017|   881|  0.5658069|       17|       35|
| V993     |   0.3547464|  0.0424604|  0.5679017|   969|  0.5658069|       17|       35|
| V670     |   0.3491045|  0.0418134|  0.5709318|   646|  0.5604856|       18|       36|
| V128     |   0.2613215|  0.0316319|  0.6230199|   105|  0.4731769|       18|       37|
| V140     |   0.2613215|  0.0316319|  0.6230199|   116|  0.4731769|       18|       37|
| V455     |   0.2613215|  0.0316319|  0.6230199|   431|  0.4731769|       18|       37|
| V511     |   0.2613215|  0.0316319|  0.6230199|   487|  0.4731769|       18|       37|
| V885     |   0.2613215|  0.0316319|  0.6230199|   861|  0.4731769|       18|       37|
| V274     |   0.2085167|  0.0254025|  0.6600751|   250|  0.4154016|       19|       38|
| V302     |   0.2085167|  0.0254025|  0.6600751|   278|  0.4154016|       19|       38|
| V552     |   0.2085167|  0.0254025|  0.6600751|   528|  0.4154016|       19|       38|
| V553     |   0.2085167|  0.0254025|  0.6600751|   529|  0.4154016|       19|       38|
| V487     |   0.2075495|  0.0252876|  0.6608057|   463|  0.4142954|       20|       39|
| V592     |   0.2075495|  0.0252876|  0.6608057|   568|  0.4142954|       20|       39|
| V663     |   0.2075495|  0.0252876|  0.6608057|   639|  0.4142954|       20|       39|
| V685     |   0.0758281|  0.0093895|  0.7900115|   661|  0.2357078|       21|       40|
| V154     |   0.0327876|  0.0040817|  0.8608133|   130|  0.1498776|       22|       41|
| V587     |   0.0249397|  0.0031078|  0.8784308|   563|  0.1296181|       22|       42|
| V195     |   0.0210763|  0.0026276|  0.8881624|   171|  0.1186007|       23|       43|
| V1026    |   0.0175782|  0.0021925|  0.8977975|  1002|  0.1078108|       23|       44|
| V137     |   0.0083076|  0.0010374|  0.9296174|   114|  0.0729822|       24|       45|

``` r
# Manhattan plot-like plot showing the effect per site colored by colinear group
ggplot(df_anova,aes(x=factor(site),y=LogP,fill=group80)) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = group_cols) +
  scale_fill_gradientn(colours = terrain.colors(27)) +
  xlab("Polymorphic sites") +
  ylab(paste("-log(P)")) +
  theme(
    panel.background = element_rect(fill = "white",colour = "white"),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid.major = element_line(),
    axis.text.x = element_text(face="bold", size=5,angle = 90, vjust = 0.5, hjust=1),
    axis.text.y = element_text(face="bold", size=10), 
    axis.title = element_text(face="bold", size=13)) +
  geom_hline(yintercept = -log(0.01),color="red",linetype='dotted') +
  geom_hline(yintercept = -log(0.05),color="blue",linetype='dotted') +
  annotate("text",x=7,y=4.7,label="P < 0.01") + 
  annotate("text",x=7,y=3.1,label="P < 0.05")
```

![](Epistasis_few_sites_files/figure-markdown_github/manplot-1.png)

Model selection, Part I
-----------------------

Build increasingly nested ANOVA models (by marginal R^2) and select best model based on LRT and AIC.

``` r
# Select representative sites per colinear group organized by marginal R^2
sorted_sites_by_R2 <- arrange(df_anova,desc(R2))$position
independent_sorted_sites_by_R2_g80 <- sorted_sites_by_R2[sorted_sites_by_R2 %in% colnames(matrix_anova_indep_1)]
independent_sorted_sites_by_R2_g99 <- sorted_sites_by_R2[sorted_sites_by_R2 %in% colnames(matrix_anova_indep_2)]

# Compute model statistics for increasingly nested linear models: Each model is a linear combination of all sites up to site (group) i (sites with an absolute correlation > 0.8)
r2_per_group <- c()
adj_r2 <- c()
group <- c()
var <- c()
aic <- c()
formulas <- list()
for(i in 1:length(independent_sorted_sites_by_R2_g80)){
  var <- c(var,independent_sorted_sites_by_R2_g80[i])
  exp_var <- sapply(var, function(x) paste0("+",x))
  formula = as.formula(paste0("effect", " ~ ", "covariate", paste0(exp_var,collapse = " ")))
  formulas[[i]] <- formula
  m <- aov(formula,data=matrix_anova_indep_1)
  r2_per_group <- c(r2_per_group,get_R2(m))
  adj_r2 <- c(adj_r2,adj_R2(m))
  aic <- c(aic,AIC(m))
  group <- c(group,i)
}

# LRT for all models
p_LRT <- c(NA)
for(i in 2:length(formulas)-1){
  pval <- anova(aov(formulas[[i]],data=matrix_anova_indep_1),aov(formulas[[i+1]],data=matrix_anova_indep_1),test="LRT")[[5]][2]
  p_LRT <- c(p_LRT,pval)
}

model <- make.unique(rep(letters, length.out = length(group)), sep='')
df_models_g80 <- data.frame(model,group,r2_per_group,adj_r2,p_LRT,aic)

# Show first 10 models:
knitr::kable(head(df_models_g80,10))
```

| model |  group|  r2\_per\_group|    adj\_r2|     p\_LRT|        aic|
|:------|------:|---------------:|----------:|----------:|----------:|
| a     |      1|       0.5995211|  0.4994014|         NA|  121.35140|
| b     |      2|       0.7752917|  0.6789881|  0.0192844|  116.99495|
| c     |      3|       0.8007926|  0.6679877|  0.3808127|  117.66991|
| d     |      4|       0.8272430|  0.6544859|  0.3816004|  118.10284|
| e     |      5|       0.8339467|  0.5848666|  0.6877945|  119.66749|
| f     |      6|       0.8508927|  0.5029757|  0.5592801|  120.48342|
| g     |      7|       0.9060100|  0.5300498|  0.2788178|  117.40713|
| h     |      8|       0.9986534|  0.9865344|  0.0000000|   72.69779|
| i     |      9|       0.9986534|  0.9865344|         NA|   72.69779|
| j     |     10|       0.9986534|  0.9865344|         NA|   72.69779|

``` r
# Compute model statistics for increasingly nested linear models: Each model is a linear combination of all sites up to site (group) i (sites with an absolute correlation > 0.8)
r2_per_group <- c()
adj_r2 <- c()
group <- c()
var <- c()
aic <- c()
formulas <- list()
for(i in 1:length(independent_sorted_sites_by_R2_g99)){
  var <- c(var,independent_sorted_sites_by_R2_g99[i])
  exp_var <- sapply(var, function(x) paste0("+",x))
  formula = as.formula(paste0("effect", " ~ ", "covariate", paste0(exp_var,collapse = " ")))
  formulas[[i]] <- formula
  m <- aov(formula,data=matrix_anova_indep_2)
  r2_per_group <- c(r2_per_group,get_R2(m))
  adj_r2 <- c(adj_r2,adj_R2(m))
  aic <- c(aic,AIC(m))
  group <- c(group,i)
}

# LRT for all models
p_LRT <- c(NA)
for(i in 2:length(formulas)-1){
  pval <- anova(aov(formulas[[i]],data=matrix_anova_indep_2),aov(formulas[[i+1]],data=matrix_anova_indep_2),test="LRT")[[5]][2]
  p_LRT <- c(p_LRT,pval)
}

model <- make.unique(rep(letters, length.out = length(group)), sep='')
df_models_g99 <- data.frame(model,group,r2_per_group,adj_r2,p_LRT,aic)

# Show first 10 models:
knitr::kable(head(df_models_g99,10))
```

| model |  group|  r2\_per\_group|    adj\_r2|     p\_LRT|       aic|
|:------|------:|---------------:|----------:|----------:|---------:|
| a     |      1|       0.6626033|  0.5782541|         NA|  119.4660|
| b     |      2|       0.6850202|  0.5500289|  0.4802979|  120.7097|
| c     |      3|       0.7138161|  0.5230268|  0.4371610|  121.6551|
| d     |      4|       0.8272430|  0.6544859|  0.0700068|  118.1028|
| e     |      5|       0.8339467|  0.5848666|  0.6877945|  119.6675|
| f     |      6|       0.8339467|  0.5848666|         NA|  119.6675|
| g     |      7|       0.8339467|  0.5848666|         NA|  119.6675|
| h     |      8|       0.8339467|  0.5848666|         NA|  119.6675|
| i     |      9|       0.8619802|  0.5399340|  0.4350357|  119.6335|
| j     |     10|       0.8619802|  0.5399340|         NA|  119.6335|

``` r
# Plot results
label1a <- "Best model: Model B\n (groups 1-2: 16 sites)"
label1b <- expression(paste("R"^2,"= 0.775"))
p1 <- ggplot(df_models_g80,aes(x=group,y=r2_per_group)) + geom_line(size=1.3) + theme_classic() +
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10), 
        axis.title = element_text(face="bold", size=13)) + 
  xlab("Group of colinear sites") + 
  ylab(expression(paste("R"^2," of nested models"))) + geom_vline(xintercept = 2,linetype="dashed",color="red") +
  ggtitle(expression(paste("Sites grouped by ", "|",rho,"|", " >0.8")))  + 
  annotate("text",x=11,y=0.75,label=label1a) + annotate("text",x=11,y=0.71,label=label1b)

label2a <- "Best model: Model D\n (groups 1-4: 16 sites)"
label2b <- expression(paste("R"^2,"= 0.827"))
p2 <- ggplot(df_models_g99,aes(x=group,y=r2_per_group)) + geom_line(size=1.3) + theme_classic() +
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10), 
        axis.title = element_text(face="bold", size=13)) + 
  xlab("Group of colinear sites") + 
  ylab(expression(paste("R"^2," of nested models"))) + geom_vline(xintercept = 4,linetype="dashed",color="red") +
  ggtitle(expression(paste("Sites grouped by ", "|",rho,"|", " >0.99"))) + 
  annotate("text",x=20,y=0.8,label=label2a) + annotate("text",x=20,y=0.76,label=label2b)

p1 + p2
```

![](Epistasis_few_sites_files/figure-markdown_github/modsel-1.png)

Model selection, Part II
------------------------

Plot the correlation between the divergence at the number of sites per group and the difference in functional effects, as a function of the cumulative number of sites per group.

``` r
##Effects of divergence at sites in group 1 through X --> Fit correlations for each nested model
correlation_per_group <- function(df,grouping,bootstrap){
  if(grouping == 80){
    groups <- unique(df$group80)
    column <- df$group80
  }
  if(grouping == 99){
    groups <- unique(df$group99)
    column <- df$group99
  }
  
  rho_cor <- c()
  number_sites <- c()
  p_val <- c()
  
  for(i in 1:length(groups)){
    group.df <- df[column %in% 1:groups[i],]
    group.sites <- as.numeric(gsub("V","",arrange(group.df,desc(R2))$position))
    # New cumulative subalignment
    group_aln <- filter_alignment_by_taxa_and_sites(taxa,group.sites,aln)
    pairwise_dists_group_sites <- get_pairwise_dist(group_aln,df_effect)
    number_sites <- c(number_sites,dim(group_aln)[2])
    m <- cor.test(df_effect$PercentDiff,pairwise_dists_group_sites[[1]],method="pearson")
    rho_cor <- c(rho_cor,m$estimate)
    
    # Compute a p-value for each correlation based on permutation
    null_cor_AAdist <- c()
    for(j in 1:bootstrap){
      perm_effects <- sample(df_effect$PercentDiff,size = length(df_effect$PercentDiff),replace = T)
      m.boot <- cor.test(pairwise_dists_group_sites[[1]],perm_effects,method="pearson")
      null_cor_AAdist <- c(null_cor_AAdist,m.boot$estimate)
    }
    p <- sum(null_cor_AAdist >= m$estimate) / bootstrap
    p_val <- c(p_val,p)
  }
  data <- data.frame(groups,number_sites,rho_cor,p_val)
  data
}

correlation_group80 <- correlation_per_group(df = df_anova,grouping = 80,bootstrap = 1000)
correlation_group99 <- correlation_per_group(df = df_anova,grouping = 99, bootstrap = 1000)


#Ploting correlations as a function of the cumulative number of sites per group
p1 <- correlation_group80 %>% dplyr::filter(.,groups %in% 1:8) %>%
  ggplot(.,aes(x=groups,y=rho_cor)) + 
  geom_point(shape=21,fill="black",aes(size=-log(p_val))) +
  #scale_fill_gradientn(colours = viridis::viridis(20)) +
  geom_line(aes(x=groups,y=rho_cor)) +
  geom_bar(aes(x=groups,y=number_sites*0.0178), stat="identity", size=.1, color="black", alpha=.4) + # 0.013 = max(rho) / max(number_sites)
  labs(x="Group",y=expression(paste(rho,": divergence at sites and effect")),size="-log(P-value)") +
  theme(
    panel.background = element_rect(fill = "white",colour = "white"),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid.major = element_line(),
    axis.text.x = element_text(face="bold", size=13),
    axis.text.y = element_text(face="bold", size=13),
    axis.title.x = element_text(face="bold", size=14),
    axis.title.y = element_text(face="bold", size=14),
    #axis.title = element_text(face="bold", size=13),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size=8),
    legend.title = element_text(size=12),
    legend.position = "top",
    legend.direction = "horizontal") +
  scale_y_continuous(sec.axis = sec_axis(~.,breaks=c(0.178,0.356,0.534,0.712,0.89), # e.g. 10*0.0178 = 0.178
                                         labels=c("10", "20", "30","40","50"),
                                         name="Cumulative number of sites")) +
  ggtitle(expression(paste("Sites grouped by ", "|",rho,"|", " >0.8")))


p2 <- correlation_group99 %>% dplyr::filter(.,groups %in% 1:12) %>%
  ggplot(.,aes(x=groups,y=rho_cor)) + 
  geom_point(shape=21,fill="black",aes(size=-log(p_val))) +
  #scale_fill_gradientn(colours = viridis::viridis(20)) +
  geom_line(aes(x=groups,y=rho_cor)) +
  geom_bar(aes(x=groups,y=number_sites*0.026), stat="identity", size=.1, color="black", alpha=.4) + # 0.013 = max(rho) / max(number_sites)
  labs(x="Group",y=expression(paste(rho,": divergence at sites and effect")),size="-log(P-value)") +
  theme(
    panel.background = element_rect(fill = "white",colour = "white"),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid.major = element_line(),
    axis.text.x = element_text(face="bold", size=13),
    axis.text.y = element_text(face="bold", size=13),
    axis.title.x = element_text(face="bold", size=14),
    axis.title.y = element_text(face="bold", size=14),
    #axis.title = element_text(face="bold", size=13),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size=8),
    legend.title = element_text(size=12),
    legend.position = "top",
    legend.direction = "horizontal") +
  scale_y_continuous(sec.axis = sec_axis(~.,breaks=c(0.26,0.52,0.78,1.04), # e.g. 10*0.0178 = 0.178
                                         labels=c("10", "20", "30","40"),
                                         name="Cumulative number of sites")) +
  ggtitle(expression(paste("Sites grouped by ", "|",rho,"|", " >0.99")))


p1 + p2
```

![](Epistasis_few_sites_files/figure-markdown_github/modesel2-1.png)

These figures show that the correlation between the divergence at the sites in the group and the difference in effects on activity (line) decreases with increasing number of sites (bars).

The effect of divergence at the 16 sites
----------------------------------------

Replot the first figure but with the divergence at the 16 sites that we identified using ANOVA.

``` r
# Compute divergence at the 16 sites 
group_sites <- as.numeric(gsub("V","",arrange(df_anova,desc(R2))[1:16,]$position))
group_aln <- filter_alignment_by_taxa_and_sites(taxa,group_sites,aln)
pairwise_dists_group_sites <- get_pairwise_dist(group_aln,df_effect)
df_effect$DistanceFew <- pairwise_dists_group_sites[[1]]

# Correlation
cor.test(df_effect$DistanceFew,df_effect$PercentDiff)
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  df_effect$DistanceFew and df_effect$PercentDiff
    ## t = 3.8392, df = 9, p-value = 0.003971
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.3566871 0.9423900
    ## sample estimates:
    ##       cor 
    ## 0.7879634

``` r
# Plot
ggplot(df_effect,aes(x=DistanceFew,y=PercentDiff),color=AAState) + 
  geom_point(size=4,shape=21,aes(fill = AAState)) +
  scale_fill_manual(values = my_pallette) +
  geom_smooth(method = lm,se = FALSE,col="grey80") +
  labs(x="# Amino acid differences between ATP1A1 orthologs \n (16 sites ANOVA model)",y="Difference in effect on two backgrounds (%)",fill="AA State") +
  theme(
    panel.background = element_rect(fill = "white",colour = "white"),
    panel.border = element_rect(linetype = 1, fill = NA),
    panel.grid.major = element_line(),
    axis.text.x = element_text(face="bold", size=10),
    axis.text.y = element_text(face="bold", size=10), 
    axis.title = element_text(face="bold", size=13),
    legend.box.background = element_blank(),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size=16),
    legend.title = element_text(size=16)) +
  annotate("text",x=5,y=180,label=expression(paste(rho, " = 0.78, P < 0.01")))
```

    ## `geom_smooth()` using formula 'y ~ x'

![](Epistasis_few_sites_files/figure-markdown_github/16siteseffect-1.png)

Permutation test
----------------

In the Model Selection, Part II section, we showed that the pattern is unlikely to be an artifact based on permutations. Here we are plotting a null distribution of correlations and comparing it with the observed correlation (making sure that the signal we're getting from the 16 sites is not an artifact).

``` r
# Repeated analyses on datasets that are permutations of the effect sizes on the constructs. 
# This would give a null distribution for the procedure and show that the amount of variation we 
# are explaining (nested ANOVA) with the sites/groups is not an artifact of the statistical procedure.

effects <- matrix_anova_indep_1$effect
cov_g1_g2 <- c("covariate","V1003","V303") # representative sites of groups 1 and 2
matrix_g1.2 <-  matrix_anova_indep_1[,colnames( matrix_anova_indep_1) %in% cov_g1_g2]

# ANOVA fitted for best model (model B, see table above)
mod_group1.2 <- aov(effect ~ covariate + V1003 + V303, data =  matrix_anova_indep_1) # 16 sites

# Perform permutation
null_r2 <- c()
for(i in 1:1000){
  perm_effects <- sample(effects,size = length(effects),replace = T) # create permutated effects (decorrelate signal)
  matrix_g1.2$perm_effects <- perm_effects
  # Do ANOVAS on permutated datasets
  m <- aov(perm_effects ~ ., data = matrix_g1.2)
  cov_i <- summary(m)[[1]][["Sum Sq"]][1]
  res_i <- tail(summary(m)[[1]][["Sum Sq"]],1)
  ss_factors <- summary(m)[[1]][["Sum Sq"]][-c(which(summary(m)[[1]][["Sum Sq"]]==cov_i),which(summary(m)[[1]][["Sum Sq"]]==res_i))]
  r2 <- sum(ss_factors)/sum(summary(m)[[1]][["Sum Sq"]])
  null_r2 <- c(null_r2,r2)
}

hist(null_r2,breaks = 30,xlab=expression(paste("R^2 on permutated effects \n(joint ANOVA model)")),main="")
abline(v = get_R2(mod_group1.2) ,col="red")
```

![](Epistasis_few_sites_files/figure-markdown_github/perm-1.png)

``` r
p <- sum(null_r2>=get_R2(mod_group1.2))/length(null_r2) # pval = 0.003

print(paste("P-value of observing an R^2 >= as observed: ", p))
```

    ## [1] "P-value of observing an R^2 >= as observed:  0.009"
