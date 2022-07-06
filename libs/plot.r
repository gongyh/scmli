#!/usr/bin/env Rscript
library(ggplot2)

args<-commandArgs(T)

#df_reads
df_stats<-read.table(file=paste0(args[1],'.stats'),col.names=c('Key','Value'))
df_reads<-df_stats[c(1,2,4,3),]
df_reads$Key<-c('All reads', 'Valid reads', 'gRNA reads', 'Unknow reads')
df_reads$Key<-factor(df_reads$Key,levels = c('All reads', 'Valid reads', 'gRNA reads', 'Unknow reads'))
names(df_reads)[names(df_reads) == 'Value'] <- 'Counts'

#df_gRNAs(only gRNAs, no unknow) var,CV
df_all<-read.table(file =paste0(args[1],'.percent'),sep = "\t", header = T)
df_gRNAs<-df_all[df_all$gene_id!='unknow',]
t<-sum(df_gRNAs$counts)
df_gRNAs$gRNAs_percent<-(df_gRNAs$counts/t)*100

Mean<-mean(df_gRNAs$gRNAs_percent)
Var<-var(df_gRNAs$gRNAs_percent)
CV<-sqrt(Var)/Mean


#plot1 reads
reads = ggplot(df_reads,aes(Key, Counts)) + 
  geom_col(width = 0.7, color = 'white', fill='lightblue') +
  geom_text(aes(label=Counts),vjust = -0.2) + xlab('') +
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(size=16),
        axis.title.y= element_text(size=16),
        axis.text.y = element_text(size=12))
ggsave('reads.png',reads)

#plot2 frequency
frequency <- ggplot(df_gRNAs,aes(x=reorder(gene_id,gRNAs_percent),y=gRNAs_percent)) + 
  theme_bw() +
  xlab('Each kind of gRNAs') +  ylab('Frequency(%)') +
  scale_y_continuous(labels = scales::comma) +
  geom_point(size=0.4) +
  geom_hline(aes(yintercept = mean(gRNAs_percent)),linetype=2,color='red') +
  geom_hline(aes(yintercept = quantile(gRNAs_percent,0.75))) +
  geom_hline(aes(yintercept = quantile(gRNAs_percent,0.25))) +
  geom_hline(aes(yintercept = median(gRNAs_percent))) +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16))
ggsave('frequency.png',frequency)

#plot3 histogram
histogram<-ggplot(df_gRNAs,aes(x=gRNAs_percent))+geom_histogram(bins=30,fill='lightblue',color='black') +
  geom_vline(aes(xintercept = mean(gRNAs_percent)),linetype=2) +
  xlab('Frequency(%)') + ylab('Counts') +
  scale_x_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16))
ggsave('histogram.png',histogram)



