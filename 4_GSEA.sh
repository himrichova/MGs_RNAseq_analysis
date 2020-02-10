# GSEA analysis

# Gene set enrichment analysis of 984 genes significantly down-regulated upon THZ531 treatment (log2FC < -4 and adj. P-value < 0.05) in comparison to dCeMM2/3/4 treatments
# Gene rankings: based on their expression change (dCeMM2/3/4 treatment vs. DMSO), using the log2FC value from differential gene expression analysis
#


for i in *.lfc2.rnk
do
echo $i

java -cp /home/himrichova/software/gsea/gsea-3.0.jar -Xmx1024m -Xms512m xtools.gsea.GseaPreranked -rpt_label ${i%.rnk}_inhouse_THZ531_lfc4_padj05_DN -gmx THZ531_lfc4_padj05_DN.gmx -rnk ${i} -set_min 20 -set_max 3500 -plot_top_x 30 -nperm 1000 -create_svgs true -scoring_scheme weighted -norm meandiv

done
