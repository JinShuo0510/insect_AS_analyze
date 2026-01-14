library(Gviz)
library(GenomicFeatures)

# 构建 TxDb 对象（需 GTF/GFF 文件）
txdb <- makeTxDbFromGFF("/data/home/jinshuo/transcript/Genome/Bombyx_mori.gtf", format="gtf")

# 基因坐标
gr_gene <- genes(txdb, filter=list(gene_id=target_gene))

# 外显子轨迹
tx_track <- GeneRegionTrack(txdb, genome="Bm", chromosome=seqnames(gr_gene),
                            transcript=tx_ids, name="Isoforms",
                            transcriptAnnotation="transcript_id",
                            showId=TRUE)

# 绘制
plotTracks(tx_track, from=start(gr_gene)-1000, to=end(gr_gene)+1000,
           main=paste(target_gene, "Transcript Structure"))
