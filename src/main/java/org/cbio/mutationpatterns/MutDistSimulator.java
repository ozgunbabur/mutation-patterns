package org.cbio.mutationpatterns;

import org.cbio.causality.util.Histogram;

import java.util.*;

/**
 * Created by BaburO on 2/28/2015.
 */
public class MutDistSimulator
{
	private Random r = new Random();

	private int find(int[] cumulative, int val)
	{
		int i = 0;
		int j = cumulative.length - 1;
		int k = j / 2;

		while (!(cumulative[k] >= val && (k == 0 || cumulative[k - 1] < val)))
		{
			if (cumulative[k] > val) j = k - 1;
			else i = k + 1;

			k = (i + j) / 2;
		}

		return k;
	}

	private int[] generateGenes(int size, int mean, double std, int min)
	{
		int[] genes = new int[size];
		for (int i = 0; i < size; i++)
		{
			genes[i] = Math.max(min, (int) ((r.nextGaussian() * std) + mean));
		}
		return genes;
	}

	private int mutateOne(int[] genes, int[] cumulative)
	{
		int val = r.nextInt(cumulative[cumulative.length - 1]);
		return find(cumulative, val);
	}

	private int[] generateCumulative(int[] genes)
	{
		int[] cumulative = new int[genes.length];
		cumulative[0] = genes[0];
		for (int i = 1; i < genes.length; i++)
		{
			cumulative[i] = genes[i] + cumulative[i - 1];
		}
		return cumulative;
	}

	private int[] mutateMany(int[] genes, int[] cumulative, int howMany)
	{
		int[] hits = new int[genes.length];
		Arrays.fill(hits, 0);
		for (int i = 0; i < howMany; i++)
		{
			hits[mutateOne(genes, cumulative)]++;
		}
		return hits;
	}

	private List<Integer>[] profile(int[] genes, int[] cumulative, int mutPerSet, int trials)
	{
		List<Integer>[] prof = new List[genes.length];
		for (int i = 0; i < prof.length; i++)
		{
			prof[i] = new ArrayList<Integer>();
		}
		for (int i = 0; i < trials; i++)
		{
			int[] hits = mutateMany(genes, cumulative, mutPerSet);
			for (int j = 0; j < hits.length; j++)
			{
				prof[j].add(hits[j]);
			}
		}

		for (List<Integer> p : prof)
		{
			Collections.sort(p);
			Collections.reverse(p);
		}

		return prof;
	}

//	private List<Double>[] generatePvals(List<Integer>[] profile)
	private void generatePvals(List<Integer>[] profile)
	{
		Histogram h = new Histogram(0.1);
		h.setBorderAtZero(true);

//		List<Double>[] pvals = new List[profile.length];
		for (int i = 0; i < profile.length; i++)
		{
//			pvals[i] = new ArrayList<Double>();

			for (int j = 0; j < profile[i].size();)
			{
				int k = j + 1;
				while (k < profile[i].size() && profile[i].get(k) == profile[i].get(j)) k++;

				double pval = k / (double) profile[i].size();
				for (int l = j; l < k; l++)
				{
					h.count(pval);
//					pvals[i].add(pval);
				}
				j = k;
			}

//			assert pvals[i].size() == profile[i].size();
		}
//		return pvals;
		h.printDensity();
	}

	private void printHistogram(List<Double>[] pvals)
	{
		Histogram h = new Histogram(0.1);
		h.setBorderAtZero(true);

		for (List<Double> list : pvals)
		{
			for (Double val : list)
			{
				h.count(val);
			}
		}

		h.printDensity();
	}

	public void run(int geneSize, int mean, double std, int min, int mutPerSample, int trials)
	{
		int[] genes = generateGenes(geneSize, mean, std, min);
		System.out.println("got genes");
		int[] cumulative = generateCumulative(genes);
		System.out.println("got cumulative");
		List<Integer>[] prof = profile(genes, cumulative, mutPerSample, trials);
		System.out.println("got profile");
		generatePvals(prof);
//		printHistogram(pvals);
	}

	public static void main(String[] args)
	{
		MutDistSimulator sim = new MutDistSimulator();
		sim.run(1000, 200, 1, 10, 3000, 10000);
	}
}
