package org.cbio.mutationpatterns;

import org.cbio.causality.data.portal.CBioPortalAccessor;
import org.cbio.causality.data.portal.CBioPortalManager;
import org.cbio.causality.data.portal.ExpDataManager;
import org.cbio.causality.data.portal.GeneticProfile;
import org.cbio.causality.util.ArrayUtil;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Ozgun Babur
 */
public class Gene
{
	String name;
	String[] mut;
	String[] cnv;
	double[] exp;

	Gene(String name, CBioPortalAccessor acc, ExpDataManager expMan)
	{
		this.name = name;
		CBioPortalManager man = acc.getManager();
		for (GeneticProfile gp : acc.getCurrentGeneticProfiles())
		{
			if (gp.getId().endsWith("mutations"))
				mut = man.getDataForGene(name, gp, acc.getCurrentCaseList());
			else if (gp.getId().endsWith("gistic"))
				cnv = man.getDataForGene(name, gp, acc.getCurrentCaseList());
		}
		exp = expMan.get(name);
	}

	public boolean[] getMutated()
	{
		boolean[] b = new boolean[mut.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = !mut[i].equals("NaN");
		}
		return b;
	}

	public boolean[] getNonMutated()
	{
		boolean[] b = new boolean[mut.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = mut[i].equals("NaN");
		}
		return b;
	}

	public boolean[] getMutatedInactivating()
	{
		boolean[] b = new boolean[mut.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = isInactivatingMutation(mut[i]);
		}
		return b;
	}

	public boolean[] getMutatedMissense()
	{
		boolean[] b = new boolean[mut.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = isMissenseMutation(mut[i]);
		}
		return b;
	}

	public boolean[] getDeleted()
	{
		boolean[] b = new boolean[cnv.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = !cnv[i].equals("NaN") && Integer.parseInt(cnv[i]) < 0;
		}
		return b;
	}

	public boolean[] getAmplified()
	{
		boolean[] b = new boolean[cnv.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = !cnv[i].equals("NaN") && Integer.parseInt(cnv[i]) > 0;
		}
		return b;
	}

	public boolean[] getCNV(int val)
	{
		boolean[] b = new boolean[cnv.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = !cnv[i].equals("NaN") && Integer.parseInt(cnv[i]) == val;
		}
		return b;
	}

	private boolean isInactivatingMutation(String val)
	{
		return val.contains("*") || val.contains("fs") || val.contains("splice");
	}

	private boolean isMissenseMutation(String val)
	{
		return !val.equals("NaN") && !isInactivatingMutation(val);
	}

	public boolean[] getEdgeExpressed(double ratio, boolean top, boolean[] consider)
	{
		assert consider == null || exp.length == consider.length;

		List<Double> list = new ArrayList<Double>(exp.length);
		for (int i = 0; i < exp.length; i++)
		{
			if (consider == null || consider[i]) list.add(exp[i]);
		}
		Collections.sort(list);
		double thr = list.get((int) (list.size() * (!top ? ratio : (1 - ratio))));

		boolean[] b = new boolean[exp.length];
		for (int i = 0; i < b.length; i++)
		{
			b[i] = !(consider != null && !consider[i]) && (top ? exp[i] >= thr : exp[i] <= thr);
		}

		return b;
	}

	public double[] cropExpression(boolean[] b)
	{
		int size = ArrayUtil.countValue(b, true);
		double[] e = new double[size];
		int k = 0;
		for (int i = 0; i < b.length; i++)
		{
			if (b[i]) e[k++] = exp[i];
		}
		return e;
	}
}
