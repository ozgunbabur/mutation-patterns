package org.cbio.mutationpatterns;

/**
 * A dataset from the portal.
 * @author Ozgun Babur
 */
public class PortalDataset implements Cloneable
{
	String name;
	String study;
	String caseList;
	String[] profile;
	boolean[] hyper;

	/**
	 * Constructor with parameters.
	 * @param study index of cancer study in portal
	 * @param caseList index of case-list in portal
	 * @param profile indices of the desired genomic profiles
	 */
	PortalDataset(String name, String study, String caseList, String[] profile)
	{
		this.name = name;
		this.study = study;
		this.caseList = caseList;
		this.profile = profile;
	}

	@Override
	protected Object clone() throws CloneNotSupportedException
	{
		PortalDataset cln = (PortalDataset) super.clone();
		if (profile != null) cln.profile = profile.clone();
		if (hyper != null) cln.hyper = hyper.clone();
		return cln;
	}
}
