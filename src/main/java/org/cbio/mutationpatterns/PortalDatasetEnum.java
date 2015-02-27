package org.cbio.mutationpatterns;

/**
 * A dataset from the portal.
 * @author Ozgun Babur
 */
public enum PortalDatasetEnum
{
	// TCGA datasets

	GBM("GBM", "gbm_tcga", "gbm_tcga_3way_complete", new String[]{"gbm_tcga_mutations", "gbm_tcga_gistic", "gbm_tcga_rna_seq_v2_mrna"}),

	OV("OV", "ov_tcga", "ov_tcga_3way_complete", new String[]{"ov_tcga_mutations", "ov_tcga_gistic", "ov_tcga_rna_seq_v2_mrna"}),

	BRCA("BRCA", "brca_tcga", "brca_tcga_3way_complete", new String[]{"brca_tcga_mutations", "brca_tcga_gistic", "brca_tcga_rna_seq_v2_mrna"}),

	COADREAD("COADREAD", "coadread_tcga_pub", "coadread_tcga_pub_3way_complete", new String[]{"coadread_tcga_pub_mutations", "coadread_tcga_pub_gistic", "coadread_tcga_pub_rna_seq_mrna"}),

	UCEC("UCEC", "ucec_tcga_pub", "ucec_tcga_pub_3way_complete", new String[]{"ucec_tcga_pub_mutations", "ucec_tcga_pub_gistic", "ucec_tcga_pub_rna_seq_v2_mrna"}),

	THCA("THCA", "thca_tcga", "thca_tcga_3way_complete", new String[]{"thca_tcga_mutations", "thca_tcga_gistic", "thca_tcga_rna_seq_v2_mrna"}),

	STAD("STAD", "stad_tcga_pub", "stad_tcga_pub_3way_complete", new String[]{"stad_tcga_pub_mutations", "stad_tcga_pub_gistic", "stad_tcga_pub_rna_seq_v2_mrna"}),

	SKCM("SKCM", "skcm_tcga", "skcm_tcga_3way_complete", new String[]{"skcm_tcga_mutations", "skcm_tcga_gistic", "skcm_tcga_rna_seq_v2_mrna"}),

	LAML("LAML", "laml_tcga", "laml_tcga_3way_complete", new String[]{"laml_tcga_mutations", "laml_tcga_gistic", "laml_tcga_rna_seq_v2_mrna"}),

	ACC("ACC", "acc_tcga", "acc_tcga_3way_complete", new String[]{"acc_tcga_mutations", "acc_tcga_gistic", "acc_tcga_rna_seq_v2_mrna"}),

	LGG("LGG", "lgg_tcga", "lgg_tcga_3way_complete", new String[]{"lgg_tcga_mutations", "lgg_tcga_gistic", "lgg_tcga_rna_seq_v2_mrna"}),

	HNSC("HNSC", "hnsc_tcga", "hnsc_tcga_3way_complete", new String[]{"hnsc_tcga_mutations", "hnsc_tcga_gistic", "hnsc_tcga_rna_seq_v2_mrna"}),

	LUAD("LUAD", "luad_tcga", "luad_tcga_3way_complete", new String[]{"luad_tcga_mutations", "luad_tcga_gistic", "luad_tcga_rna_seq_v2_mrna"}),

	LUSC("LUSC", "lusc_tcga_pub", "lusc_tcga_pub_3way_complete", new String[]{"lusc_tcga_pub_mutations", "lusc_tcga_pub_gistic", "lusc_tcga_pub_rna_seq_mrna"}),

	KIRC("KIRC", "kirc_tcga", "kirc_tcga_3way_complete", new String[]{"kirc_tcga_mutations", "kirc_tcga_gistic", "kirc_tcga_rna_seq_v2_mrna"}),

	KIRP("KIRP", "kirp_tcga", "kirp_tcga_3way_complete", new String[]{"kirp_tcga_mutations", "kirp_tcga_gistic", "kirp_tcga_rna_seq_v2_mrna"}),

	KICH("KICH", "kich_tcga", "kich_tcga_3way_complete", new String[]{"kich_tcga_mutations", "kich_tcga_gistic", "kich_tcga_rna_seq_v2_mrna"}),

	PRAD("PRAD", "prad_tcga", "prad_tcga_3way_complete", new String[]{"prad_tcga_mutations", "prad_tcga_gistic", "prad_tcga_rna_seq_v2_mrna"}),

	CESC("CESC", "cesc_tcga", "cesc_tcga_3way_complete", new String[]{"cesc_tcga_mutations", "cesc_tcga_gistic", "cesc_tcga_rna_seq_v2_mrna"}),

	LIHC("LIHC", "lihc_tcga", "lihc_tcga_3way_complete", new String[]{"lihc_tcga_mutations", "lihc_tcga_gistic", "lihc_tcga_rna_seq_v2_mrna"}),

	BLCA("BLCA", "blca_tcga", "blca_tcga_3way_complete", new String[]{"blca_tcga_mutations", "blca_tcga_gistic", "blca_tcga_rna_seq_v2_mrna"}),

	UCS("UCS", "ucs_tcga", "ucs_tcga_3way_complete", new String[]{"ucs_tcga_mutations", "ucs_tcga_gistic", "ucs_tcga_rna_seq_v2_mrna"});

	PortalDataset data;

	/**
	 * Constructor with parameters.
	 * @param study index of cancer study in portal
	 * @param caseList index of case-list in portal
	 * @param profile indices of the desired genomic profiles
	 */
	PortalDatasetEnum(String name, String study, String caseList, String[] profile)
	{
		this.data = new PortalDataset(name, study, caseList, profile);
	}

	public static PortalDataset find(String name)
	{
		for (PortalDatasetEnum data : values())
		{
			if (data.data.name.equals(name)) return data.data;
		}
		return null;
	}

	public static PortalDataset findByStudyID(String studyID)
	{
		for (PortalDatasetEnum data : values())
		{
			if (data.data.study.equals(studyID)) return data.data;
		}
		return null;
	}
}
