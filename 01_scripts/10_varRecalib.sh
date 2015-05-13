#!/bin/bash
REFERENCE="05_trinity_output/Trinity.fasta"
SNP_FOLDER="./08_callSNPs"

#Variant detection:
java -Xmx70g -jar /project/lbernatchez/drobo/users/bensuth/programs/GenomeAnalysisTK-1.0-6150-g32730b1/GenomeAnalysisTK.jar \
	-R $REFERENCE -T UnifiedGenotyper \
	-I $SNP_FOLDER/merged_realigned.bam -o $SNP_FOLDER/rawSNPS_Q4.vcf \
	-gt_mode DISCOVERY \
	-stand_call_conf 4 -stand_emit_conf 3
	
#Selecting an appropriate quality score threshold (From the Broad Institute's Wiki site):
#A common question is the confidence score threshold to use for variant detection. We recommend:
#Deep (> 10x coverage per sample) data 
#    we recommend a minimum confidence score threshold of Q30 with an emission threshold of Q10. These Q10-Q30 calls will be emitted filtered out as LowQual. 
#Shallow (< 10x coverage per sample) data 
#    because variants have by necessity lower quality with shallower coverage, we recommend a min. confidence score of Q4 and an emission threshold of Q3. 


#Annotate variants:
java -Xmx70g -jar /project/lbernatchez/drobo/users/bensuth/programs/GenomeAnalysisTK-1.0-6150-g32730b1/GenomeAnalysisTK.jar \
  -T VariantAnnotator \
  -l INFO \
  -R $REFERENCE \
  -I $SNP_FOLDER/merged_realigned.bam \
  -o $SNP_FOLDER/rawSNPS_Q4_annotated.vcf \
  -B:variant,VCF rawSNPS_Q4.vcf \
  --useAllAnnotations

#Calling InDels (needed for filtering around InDels):
java -Xmx70g -jar /project/lbernatchez/drobo/users/bensuth/programs/GenomeAnalysisTK-1.0-6150-g32730b1/GenomeAnalysisTK.jar \
	-R $REFERENCE -T UnifiedGenotyper \
	-I $SNP_FOLDER/merged_realigned.bam -o $SNP_FOLDER/InDels_Q4.vcf \
	-gt_mode DISCOVERY \
	-glm INDEL \
	-stand_call_conf 4 -stand_emit_conf 3

#Filtering around InDels:
java -Xmx70g -jar /project/lbernatchez/drobo/users/bensuth/programs/GenomeAnalysisTK-1.0-6150-g32730b1/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R $REFERENCE \
  -B:mask,VCF $SNP_FOLDER/InDels_Q4.vcf \
  -B:variant,VCF $SNP_FOLDER/rawSNPS_Q4_annotated.vcf \
  -o $SNP_FOLDER/Indel_filtered_Q4.vcf


#Additional filtering:
java -Xmx70g -jar /project/lbernatchez/drobo/users/bensuth/programs/GenomeAnalysisTK-1.0-6150-g32730b1/GenomeAnalysisTK.jar \
  -T VariantFiltration \
  -R $REFERENCE \
  -B:variant,VCF $SNP_FOLDER/Indel_filtered_Q4.vcf \
  -o $SNP_FOLDER/filtered_Q4.vcf \
  --clusterWindowSize 7 \
  --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
  --filterName "HARD_TO_VALIDATE" \
  --filterExpression "SB >= -1.0" \
  --filterName "StrandBiasFilter" \
  --filterExpression "QUAL < 10" \
  --filterName "QualFilter"
  
  
# Variant recalibrator

java -Xmx70g -jar /project/lbernatchez/drobo/users/bensuth/programs/GenomeAnalysisTK-1.0-6150-g32730b1/GenomeAnalysisTK.jar \
   -T VariantRecalibrator \
   -R $REFERENCE \
   -B:input,VCF $SNP_FOLDER/filtered_Q4.vcf \
   -B:concordantSet,VCF,known=true,training=true,truth=true,prior=10.0 $SNP_FOLDER/highQualSNPS.vcf \
   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
   -recalFile $SNP_FOLDER/VQSR.recal \
   -tranchesFile $SNP_FOLDER/VQSR.tranches \
   -rscriptFile $SNP_FOLDER/VQSR.plots.R \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
   --ignore_filter HARD_TO_VALIDATE \
   --ignore_filter LowQual
   
#Applying the recalibration
  
java -Xmx70g -jar /project/lbernatchez/drobo/users/bensuth/programs/GenomeAnalysisTK-1.0-6150-g32730b1/GenomeAnalysisTK.jar \
   -T ApplyRecalibration \
   -R $REFERENCE \
   -B:input,VCF $SNP_FOLDER/filtered_Q4.vcf \
   --ts_filter_level 99.0 \
   --ignore_filter HARD_TO_VALIDATE \
   --ignore_filter LowQual \
   -tranchesFile $SNP_FOLDER/VQSR.tranches \
   -recalFile $SNP_FOLDER/VQSR.recal \
   -o $SNP_FOLDER/recalibrated_filtered_SNPS.vcf
   
#Finally, save all the SNPS that have passed the VQSR filter into a new vcf file:

cat $SNP_FOLDER/recalibrated_filtered_SNPS.vcf | grep 'VQSLOD\|^#' | grep -v TruthSensitivityTranche > $SNP_FOLDER/VQSR_PASS_SNPS.vcf

