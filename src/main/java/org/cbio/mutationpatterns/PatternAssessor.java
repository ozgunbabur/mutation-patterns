package org.cbio.mutationpatterns;

import org.cbio.causality.util.ArrayUtil;
import org.cbio.causality.util.FDR;
import org.cbio.causality.util.FishersCombinedProbability;
import org.cbio.causality.util.Overlap;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Ozgun Babur
 */
public class PatternAssessor
{
	public static double findCNVBiasedMutations(Gene gene)
	{
		boolean[] mut = gene.getMutated();
		boolean[] oth = gene.getDeleted();

		try
		{
			return Overlap.calcCoocPval(mut, oth);
		}
		catch (Throwable t)
		{
			System.out.println();
		}
		return 0;
	}

	public static double findExpressionBiasedMutations(Gene gene)
	{
		boolean[] diploid = gene.getCNV(0);
//		boolean[] hetloss = gene.getCNV(-1);

		boolean[] mut = gene.getMutatedMissense();
		boolean[] oth1 = ArrayUtil.countValue(diploid, true) > 0 ? gene.getEdgeExpressed(0.33, true, diploid) : new boolean[diploid.length]; // over-expressed
//		boolean[] oth2 = ArrayUtil.countValue(hetloss, true) > 0 ? gene.getEdgeExpressed(0.33, true, hetloss) : new boolean[hetloss.length]; // over-expressed

//		boolean[] select = OR(diploid, hetloss);
//		boolean[] oth = OR(oth1, oth2);


//		boolean[] oth = gene.getEdgeExpressed(0.33, false, diploid); // under-expressed

//		return Overlap.calcCoocPvalOfSubset(select, mut, oth);
		return Overlap.calcCoocPvalOfSubset(diploid, mut, oth1);
	}

	private static boolean[] OR(boolean[]... b)
	{
		boolean[] r = new boolean[b[0].length];
		for (int i = 0; i < b[0].length; i++)
		{
			for (boolean[] bb : b)
			{
				if (bb[i])
				{
					r[i] = true;
					break;
				}
			}
		}
		return r;
	}

	public static void main(String[] args) throws IOException
	{
		double fdrThr = 0.1;
		for (PortalDatasetEnum dataEnum : PortalDatasetEnum.values())
		{
//			if (dataEnum != PortalDatasetEnum.LUSC) continue;

			GeneLoader loader = new GeneLoader(dataEnum.data);

			System.out.println("dataEnum = " + dataEnum);
			Map<String, Double> mutSigVals = MutSigReader.readGenesSignificance(dataEnum.name(), false);
			Map<String, Double> patternVals = new HashMap<String, Double>();

			for (String name : mutSigVals.keySet())
			{
				if (mutSigVals.get(name) == 1) continue;

				Gene gene = loader.load(name);
				if (gene != null)
				{
					patternVals.put(name, findCNVBiasedMutations(gene));
//					patternVals.put(name, findExpressionBiasedMutations(gene));
				}
			}

			Map<String, Double> combined = new HashMap<String, Double>();

			for (String name : patternVals.keySet())
			{
				combined.put(name, FishersCombinedProbability.pValue(mutSigVals.get(name), patternVals.get(name)));
			}

			List<String> selectMutSig = FDR.select(mutSigVals, null, fdrThr);
			List<String> selectPattern = FDR.select(patternVals, null, fdrThr);
			List<String> selectCombined = FDR.select(combined, null, fdrThr);

			System.out.println("selectMutSig   = " + selectMutSig);
			System.out.println("selectPattern  = " + selectPattern);
			System.out.println("selectCombined = " + selectCombined);
			selectCombined.removeAll(selectMutSig);
			System.out.println("difference     = " + selectCombined);
		}
	}
}