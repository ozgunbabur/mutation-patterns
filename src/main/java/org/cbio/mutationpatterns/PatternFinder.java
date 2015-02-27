package org.cbio.mutationpatterns;

import org.cbio.causality.data.portal.*;
import org.cbio.causality.idmapping.HGNC;
import org.cbio.causality.model.AlterationPack;
import org.cbio.causality.util.*;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class PatternFinder
{
	public static Map<String, Gene> loadAlterations(PortalDataset dataset) throws IOException
	{
		return loadAlterations(dataset, HGNC.getAllSymbols());
	}
	public static Map<String, Gene> loadAlterations(PortalDataset dataset, Set<String> symbols) throws IOException
	{
//		Set<String> symbols = new HashSet<String>(Arrays.asList("TP53"));
		CBioPortalAccessor acc = getPortalAccessor(dataset);
		ExpDataManager expMan = getExpMan(acc, dataset);
		Map<String, Gene> genes = new HashMap<String, Gene>();
//		double[] cnt = null;

		for (String symbol : symbols)
		{
			AlterationPack pack = acc.getAlterations(symbol);

			if (pack != null)
			{
//				if (cnt == null) cnt = new double[pack.getSize()];

//				Change[] changes = pack.get(Alteration.MUTATION);
//				for (int i = 0; i < changes.length; i++)
//				{
//					if (changes[i].isAltered()) cnt[i]++;
//				}
				genes.put(symbol, new Gene(symbol, acc, expMan));
			}
		}

//		boolean[] hyper = Summary.markOutliers(cnt, true);

		return genes;
	}

	public static Set<String> getBorderSymbols(PortalDataset data,
		double qValLowerLimit, double qValUpperLimit)
	{
		String code = data.name;
		if (code.contains("_")) code = code.substring(0, code.indexOf("_"));

		Set<String> genes = BroadAccessor.getMutsigGenes(code, qValUpperLimit, true);
		genes.removeAll(BroadAccessor.getMutsigGenes(code, qValLowerLimit, true));
		return genes;
	}

	public static CBioPortalAccessor getPortalAccessor(PortalDataset data) throws IOException
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
	
	public static ExpDataManager getExpMan(CBioPortalAccessor acc, PortalDataset dataset) throws IOException
	{
		GeneticProfile profile = findProfile(acc.getGeneticProfilesForCurrentStudy(), dataset.profile[2]);
		ExpDataManager eman = new ExpDataManager(profile, acc.getCurrentCaseList());
		eman.setTakeLog(true);
		return eman;
	}

	private static CancerStudy findCancerStudy(List<CancerStudy> list, String id)
	{
		for (CancerStudy study : list)
		{
			if (study.getStudyId().equals(id)) return study;
		}
		return null;
	}

	private static CaseList findCaseList(List<CaseList> list, String id)
	{
		for (CaseList cl : list)
		{
			if (cl.getId().equals(id)) return cl;
		}
		return null;
	}

	private static GeneticProfile findProfile(List<GeneticProfile> list, String id)
	{
		for (GeneticProfile profile : list)
		{
			if (profile.getId().equals(id)) return profile;
		}
		return null;
	}

	public static Map<String, Double>[] findCNVBiasedMutations(Map<String, Gene> genes)
	{
		Map<String, Double> pvals = new HashMap<String, Double>();
		Map<String, Double> limits = new HashMap<String, Double>();
		int testedCnt = 0;
		for (String name : genes.keySet())
		{
			Gene gene = genes.get(name);

			boolean[] mut = gene.getMutatedInactivating();
//			boolean[] mut = gene.getMutatedMissense();
//			boolean[] mut = gene.getMutated();

			boolean[] oth = gene.getDeleted();
//			boolean[] oth = gene.getAmplified();

			int mutCnt = ArrayUtil.countValue(mut, true);
			int othCnt = ArrayUtil.countValue(oth, true);

			if (mutCnt < 5 || othCnt < 5) continue;
			if (Summary.mean(gene.exp) < 4) continue;

			testedCnt++;

			double pv = Overlap.calcCoocPvalOfSubset(mut, oth);
			double lim = Overlap.calcCoocPval(mut.length, Math.min(mutCnt, othCnt), mutCnt, othCnt);

			pvals.put(name, pv);
			limits.put(name, lim);
		}
		List<String> select = FDR.select(pvals, limits, 0.05);
		if (!select.isEmpty())
		{
			Map<String, Double> qVals = FDR.getQVals(pvals, limits);
			for (String gene : select)
			{
				System.out.println(gene + "\t" +
					FormatUtil.roundToSignificantDigits(pvals.get(gene), 2) + "\t" +
					FormatUtil.roundToSignificantDigits(qVals.get(gene), 2));
			}
		}

		System.out.println("testedCnt = " + testedCnt);
		return new Map[]{pvals, limits};
	}

	public static Map<String, Double>[] findExpressionBiasedMutations(Map<String, Gene> genes)
	{
		Map<String, Double> pvals = new HashMap<String, Double>();
		Map<String, Double> limits = new HashMap<String, Double>();
		int testedCnt = 0;
		for (String name : genes.keySet())
		{
			Gene gene = genes.get(name);
			boolean[] diploid = gene.getDiploid();
			int n = ArrayUtil.countValue(diploid, true);
			if (n < 50) continue;

//			boolean[] mut = gene.getMutatedInactivating();
			boolean[] mut = gene.getMutatedMissense();
//			boolean[] mut = gene.getMutated();

			boolean[] oth = gene.getEdgeExpressed(0.33, true, diploid); // over-expressed
//			boolean[] oth = gene.getEdgeExpressed(0.33, false, diploid); // under-expressed

			int mutCnt = ArrayUtil.countValue(mut, diploid, true);
			int othCnt = ArrayUtil.countValue(oth, diploid, true);

			if (mutCnt < 5 || othCnt < 5) continue;
			if (Summary.mean(gene.exp) < 4) continue;

			testedCnt++;

			double pv = Overlap.calcCoocPvalOfSubset(diploid, mut, oth);
			double lim = Overlap.calcCoocPval(n, Math.min(mutCnt, othCnt), mutCnt, othCnt);

			pvals.put(name, pv);
			limits.put(name, lim);
		}
		List<String> select = FDR.select(pvals, limits, 0.05);
		if (!select.isEmpty())
		{
			Map<String, Double> qVals = FDR.getQVals(pvals, limits);
			for (String gene : select)
			{
				System.out.println(gene + "\t" +
					FormatUtil.roundToSignificantDigits(pvals.get(gene), 2) + "\t" +
					FormatUtil.roundToSignificantDigits(qVals.get(gene), 2));
			}
		}

		System.out.println("testedCnt = " + testedCnt);
		return new Map[]{pvals, limits};
	}

	private static void AND(boolean[] b1, boolean[] b2)
	{
		for (int i = 0; i < b1.length; i++)
		{
			if (!b2[i]) b1[i] = false;
		}
	}

	public static Map<String, Double>[] findBiasedGenes2(Map<String, Gene> genes)
	{
		Map<String, Double> pvals = new HashMap<String, Double>();
		Map<String, Double> changes = new HashMap<String, Double>();
		Map<String, Double> limits = new HashMap<String, Double>();
		int testedCnt = 0;
		for (String name : genes.keySet())
		{
			Gene gene = genes.get(name);

//			boolean[] mut = gene.getMutatedInactivating();
			boolean[] mut = gene.getMutatedMissense();
			boolean[] nmut = gene.getNonMutated();


			int mutCnt = ArrayUtil.countValue(mut, true);

			if (mutCnt / (double) mut.length < 0.05) continue;

			testedCnt++;

			double[] mExp = gene.cropExpression(mut);
			double[] nExp = gene.cropExpression(nmut);

			double pv = StudentsT.getPValOfMeanDifference(nExp, mExp) / 2;
			double ch = Summary.calcChangeOfMean(nExp, mExp);

			if (ch < 0) continue;

			pvals.put(name, pv);
			limits.put(name, 0D);
			changes.put(name, ch);
		}
		List<String> select = FDR.select(pvals, limits, 0.05);
		if (!select.isEmpty())
		{
			Map<String, Double> qVals = FDR.getQVals(pvals, limits);
			for (String gene : select)
			{
				System.out.println(gene + "\t" +
					FormatUtil.roundToSignificantDigits(pvals.get(gene), 2) + "\t" +
					FormatUtil.roundToSignificantDigits(qVals.get(gene), 2) + "\t" +
					(changes.get(gene) > 0 ? "+" : "-"));
			}
		}

		System.out.println("testedCnt = " + testedCnt);
		return new Map[]{pvals, limits};
	}

	private static int[] countMutDistToOther(boolean[] mut, boolean[] oth)
	{
		int t = 0;
		int f = 0;

		for (int i = 0; i < mut.length; i++)
		{
			if (mut[i])
			{
				if (oth[i]) t++; else f++;
			}
		}
		return new int[]{t, f};
	}

	public static void main(String[] args) throws IOException
	{
//		plotMutVersusExp();

		Map<String, List<Double>> vals = new HashMap<String, List<Double>>();
		Map<String, List<Double>> lims = new HashMap<String, List<Double>>();
		Map<String, List<String>> assoc = new HashMap<String, List<String>>();

		for (PortalDatasetEnum dataset : PortalDatasetEnum.values())
		{
//			if (!dataset.data.name.startsWith("HNSC")) continue;

			if (dataset.data.name.equals("simulated") || dataset.data.name.equals("sarcoma") ||
				dataset.data.name.endsWith("-pub")) continue;

			System.out.println("dataset = " + dataset.name());
//			Map<String, Gene> genes = loadAlterations(dataset.data, Collections.singleton("DNAH3"));
			Map<String, Gene> genes = loadAlterations(dataset.data);
//			Map<String, Gene> genes = loadAlterations(dataset.data, getBorderSymbols(dataset.data, 0.05, 0.99));
			Map<String, Double>[] pvs = findExpressionBiasedMutations(genes);
//			Map<String, Double>[] pvs = findBiasedGenes2(genes);
			for (String gene : pvs[0].keySet())
			{
				if (!vals.containsKey(gene)) vals.put(gene, new ArrayList<Double>());
				Double v = pvs[0].get(gene);
				vals.get(gene).add(v);
				if (v < 1)
				{
					if (!assoc.containsKey(gene)) assoc.put(gene, new ArrayList<String>());
					assoc.get(gene).add("[" + dataset.name() + " " + FormatUtil.roundToSignificantDigits(v, 2) + "]");
				}
			}
			for (String gene : pvs[1].keySet())
			{
				if (!lims.containsKey(gene)) lims.put(gene, new ArrayList<Double>());
				lims.get(gene).add(pvs[1].get(gene));
			}
			System.out.println();
		}

		System.out.println("Cumulative:");
		Map<String, Double> pvals = new HashMap<String, Double>();
		Map<String, Double> limit = new HashMap<String, Double>();
		for (String gene : vals.keySet())
		{
			pvals.put(gene, FishersCombinedProbability.pValue(ArrayUtil.toArray(vals.get(gene), 0D)));
			limit.put(gene, FishersCombinedProbability.pValue(ArrayUtil.toArray(lims.get(gene), 0D)));
		}

		double fdr = FDR.decideBestFDR_BH(pvals, limit);
		System.out.println("fdr = " + fdr);
		List<String> select = FDR.select(pvals, limit, fdr);
		System.out.println("select = " + select);
		if (!select.isEmpty())
		{
			Map<String, Double> qVals = FDR.getQVals(pvals, limit);
			for (String gene : select)
			{
				System.out.println(gene + "\t" +
					FormatUtil.roundToSignificantDigits(pvals.get(gene), 2) + "\t" +
					FormatUtil.roundToSignificantDigits(qVals.get(gene), 2) + "\t" + assoc.get(gene));
			}
		}
	}

	private static void plotMutVersusExp() throws IOException
	{
		Map<String, Gene> genes = loadAlterations(PortalDatasetEnum.GBM.data,
			new HashSet<String>(Arrays.asList("TP53", "RHOC")));

		Gene mutGene = genes.get("TP53");
		Gene expGene = genes.get("RHOC");

		boolean[] mut = mutGene.getMutated();
		int mutCnt = ArrayUtil.countValue(mut, true);

		double[] exp = expGene.exp;
		double[] exp1 = new double[mut.length - mutCnt];
		double[] exp2 = new double[mutCnt];

		int a = 0;
		int b = 0;
		for (int i = 0; i < mut.length; i++)
		{
			if (mut[i]) exp2[b++] = exp[i];
			else exp1[a++] = exp[i];
		}

		double change = Summary.calcChangeOfMean(exp1, exp2);
		double pval = StudentsT.getPValOfMeanDifference(exp1, exp2);

		System.out.println("change = " + change);
		System.out.println("pval = " + pval);

		System.exit(0);
	}

}
