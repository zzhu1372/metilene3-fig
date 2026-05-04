library(EnrichedHeatmap)
library(hash)

library(RColorBrewer)
library(circlize)  # for colorRamp2

# Use 11-color Spectral palette and reverse it (blue low, red high)
spectral_colors <- rev(brewer.pal(11, "Spectral"))

# Define color scale, assuming data range from -2 to 2
col_fun <- colorRamp2(
  breaks = seq(0, 1, length.out = 11), 
  colors = spectral_colors
)

a = read.csv('./PDAC/motif/groupmean.met.csv')
a$end = a$pos+1
b = makeGRangesFromDataFrame(a,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="pos",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)

gene = 'NFATC1'
pdf(file=paste0('./PDAC/motif/',gene,'.allhighDMRs.smooth.pdf'), width = 20, height = 4)
x <- hash()
for (i in c('mean_PDAC','mean_PanIN','mean_Acinar','mean_Duct'))
{
    c = read.table(paste0('./PDAC/motif/',gene,'.find.allhighDMRs.emap.tsv'))
    d = makeGRangesFromDataFrame(c,
                             keep.extra.columns=FALSE,
                             ignore.strand=FALSE,
                             seqinfo=NULL,
                             seqnames.field=c("V1"),
                             start.field="V2",
                             end.field=c("V3"),
                             strand.field="strand",
                             starts.in.df.are.0based=FALSE)
    
    x[[i]] = normalizeToMatrix(b, d, value_column = i, extend = 3000, w = 100, smooth = TRUE)
}
EnrichedHeatmap(x[['mean_PDAC']], col = col_fun, name = 'mean_PDAC', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )+
EnrichedHeatmap(x[['mean_PanIN']], col = col_fun, name = 'mean_PanIN', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )+
EnrichedHeatmap(x[['mean_Acinar']], col = col_fun, name = 'mean_Acinar', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )+
EnrichedHeatmap(x[['mean_Duct']], col = col_fun, name = 'mean_Duct', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )
dev.off()

gene = 'NFKB2'
pdf(file=paste0('./PDAC/motif/',gene,'.allhighDMRs.smooth.pdf'), width = 20, height = 4)
x <- hash()
for (i in c('mean_PDAC','mean_PanIN','mean_Acinar','mean_Duct'))
{
    c = read.table(paste0('./PDAC/motif/',gene,'.find.allhighDMRs.emap.tsv'))
    d = makeGRangesFromDataFrame(c,
                             keep.extra.columns=FALSE,
                             ignore.strand=FALSE,
                             seqinfo=NULL,
                             seqnames.field=c("V1"),
                             start.field="V2",
                             end.field=c("V3"),
                             strand.field="strand",
                             starts.in.df.are.0based=FALSE)
    
    x[[i]] = normalizeToMatrix(b, d, value_column = i, extend = 3000, w = 100, smooth = TRUE)
}
EnrichedHeatmap(x[['mean_PDAC']], col = col_fun, name = 'mean_PDAC', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )+
EnrichedHeatmap(x[['mean_PanIN']], col = col_fun, name = 'mean_PanIN', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )+
EnrichedHeatmap(x[['mean_Acinar']], col = col_fun, name = 'mean_Acinar', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )+
EnrichedHeatmap(x[['mean_Duct']], col = col_fun, name = 'mean_Duct', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )
dev.off()

gene = 'NFKB2_NFATC1'
pdf(file=paste0('./PDAC/motif/',gene,'.allhighDMRs.smooth.pdf'), width = 20, height = 4)
x <- hash()
for (i in c('mean_PDAC','mean_PanIN','mean_Acinar','mean_Duct'))
{
    c = read.table(paste0('./PDAC/motif/',gene,'.find.allhighDMRs.emap.tsv'))
    d = makeGRangesFromDataFrame(c,
                             keep.extra.columns=FALSE,
                             ignore.strand=FALSE,
                             seqinfo=NULL,
                             seqnames.field=c("V1"),
                             start.field="V2",
                             end.field=c("V3"),
                             strand.field="strand",
                             starts.in.df.are.0based=FALSE)
    
    x[[i]] = normalizeToMatrix(b, d, value_column = i, extend = 3000, w = 100, smooth = TRUE)
}
EnrichedHeatmap(x[['mean_PDAC']], col = col_fun, name = 'mean_PDAC', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )+
EnrichedHeatmap(x[['mean_PanIN']], col = col_fun, name = 'mean_PanIN', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )+
EnrichedHeatmap(x[['mean_Acinar']], col = col_fun, name = 'mean_Acinar', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )+
EnrichedHeatmap(x[['mean_Duct']], col = col_fun, name = 'mean_Duct', column_title = gene, top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim = c(0, 1)), height = unit(3, "cm")), )
dev.off()

print('Done.')