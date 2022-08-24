#!/usr/bin/env Rscript
library(ggplot2)

args<-commandArgs(T)

#df_reads
df_stats<-read.table(file=paste0(args[1],'.stats'),col.names=c('Key','Value'))
df_reads<-df_stats[c(2,3,5,4),]
df_reads$Key<-c('All reads', 'Valid reads', 'gRNA reads', 'Unknow reads')
df_reads$Key<-factor(df_reads$Key,levels = c('All reads', 'Valid reads', 'gRNA reads', 'Unknow reads'))
names(df_reads)[names(df_reads) == 'Value'] <- 'Counts'

#df_gRNAs(only gRNAs, no unknow)
df_all<-read.table(file =paste0(args[1],'.percentage'),sep = "\t", header = T)
df_gRNAs<-df_all[df_all$gene_id!='unknow',]
df_gRNAs_exist<-df_all[(df_all$gene_id!='unknow')&(df_all$counts!=0),]
t<-sum(df_gRNAs$counts)
df_gRNAs$gRNAs_percentage<-(df_gRNAs$counts/t)*100
t_exist<-sum(df_gRNAs_exist$counts)
df_gRNAs_exist$gRNAs_percentage<-(df_gRNAs_exist$counts/t_exist)*100

#plot1 reads
reads = ggplot(df_reads,aes(Key, Counts)) + 
  geom_col(width = 0.7, color = 'white', fill='lightblue') +
  geom_text(aes(label=Counts),vjust = -0.2) + xlab('') +
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(size=16),
        axis.title.y= element_text(size=16),
        axis.text.y = element_text(size=12))
ggsave('reads.pdf',reads)

#plot2(1/2) frequency
frequency <- ggplot(df_gRNAs,aes(x=reorder(gene_id,gRNAs_percentage),y=gRNAs_percentage)) + 
  theme_bw() +
  ylab('Frequency(%)') +
  scale_y_continuous(labels = scales::comma) +
  geom_point(size=0.4) +
  geom_hline(aes(yintercept = mean(gRNAs_percentage)),linetype=2,color='red') +
  geom_hline(aes(yintercept = quantile(gRNAs_percentage,0.75))) +
  geom_hline(aes(yintercept = quantile(gRNAs_percentage,0.25))) +
  geom_hline(aes(yintercept = median(gRNAs_percentage))) +
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16))
frequency <- frequency + aes(x=seq(1,length(gRNAs_percentage)))+scale_x_continuous()+xlab('Each kind of gRNAs')

#add text
x_range <- layer_scales(frequency)$x$range$range
y_range <- layer_scales(frequency)$y$range$range
Mean<-mean(df_gRNAs$gRNAs_percentage)
Var<-var(df_gRNAs$gRNAs_percentage)
CV<-sqrt(Var)/Mean
frequency<-frequency + geom_text(x=x_range[2]*0.8,
    y=y_range[2]*0.92,
    aes(label=paste0('Var = ',Var,'\nCV = ',CV)))
ggsave('frequency.png',frequency)

#plot2(2/2) frequency_exist
frequency_exist <- ggplot(df_gRNAs_exist,aes(x=reorder(gene_id,gRNAs_percentage),y=gRNAs_percentage)) +
  theme_bw() +
  ylab('Frequency(%)') +
  scale_y_continuous(labels = scales::comma) +
  geom_point(size=0.4) +
  geom_hline(aes(yintercept = mean(gRNAs_percentage)),linetype=2,color='red') +
  geom_hline(aes(yintercept = quantile(gRNAs_percentage,0.75))) +
  geom_hline(aes(yintercept = quantile(gRNAs_percentage,0.25))) +
  geom_hline(aes(yintercept = median(gRNAs_percentage))) +
  theme(axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16))
frequency_exist <- frequency_exist + aes(x=seq(1,length(gRNAs_percentage)))+scale_x_continuous()+xlab('Each kind of gRNAs')

#add text
x_range <- layer_scales(frequency_exist)$x$range$range
y_range <- layer_scales(frequency_exist)$y$range$range
Mean<-mean(df_gRNAs_exist$gRNAs_percentage)
Var<-var(df_gRNAs_exist$gRNAs_percentage)
CV<-sqrt(Var)/Mean
frequency_exist<-frequency_exist + geom_text(x=x_range[2]*0.8,
    y=y_range[2]*0.92,
    aes(label=paste0('Var = ',Var,'\nCV = ',CV)))
ggsave('frequency_exist.png',frequency_exist)


#plot3(1/2) histogram
histogram<-ggplot(df_gRNAs,aes(x=gRNAs_percentage))+geom_histogram(bins=30,fill='lightblue',color='black') +
  geom_vline(aes(xintercept = mean(gRNAs_percentage)),linetype=2) +
  xlab('Frequency(%)') + ylab('Counts') +
  scale_x_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16))
ggsave('histogram.pdf',histogram)

#plot3(2/2) histogram_exist
histogram_exist<-ggplot(df_gRNAs_exist,aes(x=gRNAs_percentage))+geom_histogram(bins=30,fill='lightblue',color='black') +
  geom_vline(aes(xintercept = mean(gRNAs_percentage)),linetype=2) +
  xlab('Frequency(%)') + ylab('Counts') +
  scale_x_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x= element_text(size=16),
        axis.title.y= element_text(size=16))
ggsave('histogram_exist.pdf',histogram_exist)

#plot4 cumulative_unknow_percentage
cumulative_unknow_percentage<-ggplot(df_all,aes(x=seq(1,length(percentage)),y=cumulative_unknow_percentage)) +
  geom_point(size=0.2) + xlab("Kinds of sequences") + ylab("Cumulative unknow percentage") +
  theme_bw()
ggsave('cumulative_unknow_percentage.pdf',cumulative_unknow_percentage)
