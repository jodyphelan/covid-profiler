#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}


library(ggpubr)
library(tmap)
library(svglite)
data("World")
theme_set(theme_pubr())

filename<-args[1]
forward_seq<-args[2]
reverse_seq<-args[3]
probe_Seq<-args[4]

# filename<-"~/projects/covid/analysis/primers//China_N.primer.csv"
# forward_seq<-"GGGGAACTTCTCCTGCTAGAAT"
# reverse_seq<-"CAGACATTTTGCTCTCAAGCTG"
# probe_Seq<-"TTGCTGCTGCTTGACAGATT"


csv<-read.csv(filename,stringsAsFactors = F)
csv<-subset(csv, !is.na(forward_primer_mismatches))
meta<-read.csv("~/projects/covid/scraper/round3/gisaid.csv",stringsAsFactors = F)
meta$country[which(meta$country=="USA")]<-"United States"
meta$country[which(meta$country=="England")]<-"United Kingdom"
csv$country<-as.character(meta$country[match(as.character(csv$id),meta$seqname)])


plot_mismatch_by_position<-function(column,refseq){
  z<-t(sapply(csv[,column],function(x){
    as.numeric(!apply(matrix(c(strsplit(x,"")[[1]],strsplit(refseq,"")[[1]]),byrow=T,nrow=2),2,function(z){
      match_nucleotide(z[1],z[2])
    }))
  }))


  df<-as.data.frame(apply(z,2,mean))
  df$position<-1:nrow(df)

  colnames(df)<-c("freq","position")

  ggplot(data=df, aes(x=position, y=freq)) +
    geom_line()+
    geom_point()+
    labs(y="Freq",x="Position") +
    theme(text = element_text(size=7))
}


match_nucleotide<-function(x,y){
  switch(x,
         A= ifelse(y %in% c("A","R","W","M","D","H","V","N"), TRUE, FALSE),
         C= ifelse(y %in% c("C","Y","S","M","B","H","V","N"), TRUE, FALSE),
         G= ifelse(y %in% c("G","R","S","K","B","D","V","N"), TRUE, FALSE),
         T= ifelse(y %in% c("T","Y","W","K","B","D","H","N"), TRUE, FALSE),
         R= ifelse(y %in% c("A","G"), TRUE, FALSE),
         Y= ifelse(y %in% c("C","T"), TRUE, FALSE),
         S= ifelse(y %in% c("G","C"), TRUE, FALSE),
         W= ifelse(y %in% c("A","T"), TRUE, FALSE),
         K= ifelse(y %in% c("G","T"), TRUE, FALSE),
         M= ifelse(y %in% c("A","C"), TRUE, FALSE),
         B= ifelse(y %in% c("C","G","T"), TRUE, FALSE),
         D= ifelse(y %in% c("A","G","T"), TRUE, FALSE),
         H= ifelse(y %in% c("A","C","T"), TRUE, FALSE),
         V= ifelse(y %in% c("A","C","G"), TRUE, FALSE),
         N= ifelse(y %in% c("A","C","G","T"), TRUE, FALSE)
  )
}



get_country_values<-function(col,mismatch_threshold=1){
  sapply(names(table(csv$country)),function(d){
    country_subset<-subset(csv,country==d)
    mean(country_subset[,col]>=mismatch_threshold)
  })
}


# col<-"forward_primer_mismatches"
# mismatch_threshold<-1

mismatch_barplot<-function(col){
  tmp<-as.data.frame(table(factor(csv[,col],levels=seq(0,5))))
  ggbarplot(tmp, x = "Var1", y = "Freq",
            fill = "steelblue",               # change fill color by cyl
            color = "black",
            xlab="# mismatches"
  ) +
    theme(text = element_text(size=7))
}

plot_mismatch_column<-function(col,mismatch_threshold=1){
  fp<-get_country_values(col, mismatch_threshold)

  World$Mismatch_present<-fp[match(World$name,names(fp))]

  ggplot(World) +
    geom_sf(fill = "lightgrey",lwd = 0) +
    geom_sf(data = subset(World,!is.na(Mismatch_present)),aes(fill = Mismatch_present ),lwd=0) +
    scale_fill_viridis_c(alpha=.5)  +
    guides(fill = guide_colourbar(title=NULL,label.theme = element_text(size=7)) )
}



m1<-plot_mismatch_column("forward_primer_mismatches",1)
b1<-mismatch_barplot("forward_primer_mismatches")
m2<-plot_mismatch_column("probe_mismatches",1)
b2<-mismatch_barplot("probe_mismatches")
m3<-plot_mismatch_column("reverse_primer_mismatches",1)
b3<-mismatch_barplot("reverse_primer_mismatches")



l1<-plot_mismatch_by_position("forward_primer_seq",forward_seq)
l2<-plot_mismatch_by_position("probe_sequence",probe_Seq)
l3<-plot_mismatch_by_position("reverse_primer_seq",reverse_seq)


#pdf(gsub(".csv",".pdf",filename))
svg(gsub(".csv",".svg",filename))
ggarrange(
  ggarrange(m1, ggarrange(b1,l1,nrow=2,labels = c("B","C")),ncol = 2),
  ggarrange(m2, ggarrange(b2,l2,nrow=2,labels = c("E","F")),ncol = 2),
  ggarrange(m3, ggarrange(b3,l3,nrow=2,labels = c("H","I")),ncol = 2),
  nrow=3,labels = c("A","D","G")
)
dev.off()
