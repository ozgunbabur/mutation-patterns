package org.cbio.mutationpatterns;

import org.cbio.causality.data.portal.*;
import org.cbio.causality.model.AlterationPack;

import java.io.IOException;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class GeneLoader
{
	private CBioPortalAccessor acc;
	private ExpDataManager expMan;

	public GeneLoader(PortalDataset dataset) throws IOException
	{
		acc = getPortalAccessor(dataset);
		expMan = getExpMan(acc, dataset);
	}

	public Map<String, Gene> load(Set<String> symbols) throws IOException
	{
		Map<String, Gene> genes = new HashMap<String, Gene>();

		for (String symbol : symbols)
		{
			Gene gene = load(symbol);
			if (gene != null) genes.put(symbol, gene);
		}

		return genes;
	}

	public Gene load(String symbol) throws IOException
	{
		AlterationPack pack = acc.getAlterations(symbol);

		if (pack != null)
		{
			return new Gene(symbol, acc, expMan);
		}

		return null;
	}

	public CBioPortalAccessor getPortalAccessor(PortalDataset data) throws IOException
	{
		CBioPortalAccessor acc = new CBioPortalAccessor();
		acc.setUseCacheOnly(true);
		CBioPortalOptions opts = new CBioPortalOptions();
		opts.put(CBioPortalOptions.PORTAL_OPTIONS.CNA_LOWER_THRESHOLD, -1D);
		opts.put(CBioPortalOptions.PORTAL_OPTIONS.CNA_UPPER_THRESHOLD, 1D);
		acc.setOptions(opts);

		CancerStudy cancerStudy = findCancerStudy(acc.getCancerStudies(), data.study);
		acc.setCurrentCancerStudy(cancerStudy);

		List<GeneticProfile> geneticProfilesForCurrentStudy =
			acc.getGeneticProfilesForCurrentStudy();
		List<GeneticProfile> gp = new ArrayList<GeneticProfile>();
		for (String prof : data.profile)
		{
			if (prof.endsWith("gistic") || prof.endsWith("mutations"))
			{
				gp.add(findProfile(geneticProfilesForCurrentStudy, prof));
			}
		}
		acc.setCurrentGeneticProfiles(gp);

		List<CaseList> caseLists = acc.getCaseListsForCurrentStudy();
		acc.setCurrentCaseList(findCaseList(caseLists, data.caseList));
		return acc;
	}

	public ExpDataManager getExpMan(CBioPortalAccessor acc, PortalDataset dataset) throws IOException
	{
		GeneticProfile profile = findProfile(acc.getGeneticProfilesForCurrentStudy(), dataset.profile[2]);
		ExpDataManager eman = new ExpDataManager(profile, acc.getCurrentCaseList());
		eman.setTakeLog(true);
		return eman;
	}

	private CancerStudy findCancerStudy(List<CancerStudy> list, String id)
	{
		for (CancerStudy study : list)
		{
			if (study.getStudyId().equals(id)) return study;
		}
		return null;
	}

	private CaseList findCaseList(List<CaseList> list, String id)
	{
		for (CaseList cl : list)
		{
			if (cl.getId().equals(id)) return cl;
		}
		return null;
	}

	private GeneticProfile findProfile(List<GeneticProfile> list, String id)
	{
		for (GeneticProfile profile : list)
		{
			if (profile.getId().equals(id)) return profile;
		}
		return null;
	}

}
