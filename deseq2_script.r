library(DESeq2);

counts = read.csv("<#COUNT_CSV>", row.names = 1);

colData = data.frame(condition = c("untreated", "treated"));
rownames(colData) = colnames(counts);

dds = DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition);

dds = DESeq(dds);
res = results(dds);

write.csv(as.data.frame(res), file="deseq_output.csv");
