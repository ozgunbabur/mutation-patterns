package org.cbio.mutationpatterns;

import org.cbio.causality.data.portal.*;
import org.cbio.causality.idmapping.HGNC;
import org.cbio.causality.util.*;

import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.util.*;

/**
 * @author Ozgun Babur
 */
public class MutSigReader
{
	private static final String BROAD_DIR = "broad-data/";
	private static final String CACHED_STUDIES_FILE = "studies.txt";
	private static String cacheDir;
	private static String broadDataURL = "http://gdac.broadinstitute.org/runs/analyses__latest/";
	private static List<String> studyCodes;
	private static final String MUTSIG_ANALYSIS_SUBSTR = "MutSigNozzleReport2CV.Level_4";
	public static final String[] mutsigPartialFileNames = new String[]{".sig_genes.", "sig_genes."};

	public static void setCacheDir(String dir)
	{
		cacheDir = dir;
	}

	public static void setBroadDataURL(String url)
	{
		broadDataURL = url;
		studyCodes = null;
	}

	public static String getBroadDataURL()
	{
		if (broadDataURL.endsWith("latest/"))
		{
			try
			{
				URL url = new URL(broadDataURL);
				URLConnection con = url.openConnection();
				BufferedReader reader = new BufferedReader(
					new InputStreamReader(con.getInputStream()));

				for (String line = reader.readLine(); line != null; line = reader.readLine())
				{
					if (line.startsWith("<h3>") && line.endsWith("analyses  Run</h3>"))
					{
						String date = line.substring(line.indexOf(">") + 1, line.indexOf(" "));
						System.out.println("date = " + date);
						broadDataURL = broadDataURL.substring(0, broadDataURL.lastIndexOf("l")) +
							date + "/";

						break;
					}
				}
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
		}
		return broadDataURL;
	}

	public static List<String> getStudyCodes()
	{
		if (studyCodes == null)
		{
			studyCodes = readStudiesFromCache();

			if (studyCodes == null)
			{
				studyCodes = new ArrayList<String>(30);
				try
				{
					URL url = new URL(getBroadDataURL() + "ingested_data.tsv");

					URLConnection con = url.openConnection();

					BufferedReader reader = new BufferedReader(
						new InputStreamReader(con.getInputStream()));

					for (String line = reader.readLine(); line != null; line = reader.readLine())
					{
						if (line.isEmpty() || line.startsWith("#")
							|| line.startsWith("Cohort") || line.startsWith("Totals")) continue;

						String study = line.substring(0, line.indexOf("\t"));
						if (mutSigAvailable(study)) studyCodes.add(study);
					}
					reader.close();

					// Keep only the ones that are available in cBioPortal

					Set<String> available = new HashSet<String>();
					CBioPortalAccessor acc = new CBioPortalAccessor();
					for (CancerStudy cancerStudy : acc.getCancerStudies())
					{
						if (cancerStudy.getStudyId().endsWith("_tcga"))
						{
							available.add(cancerStudy.getStudyId().substring(0,
								cancerStudy.getStudyId().indexOf("_")).toUpperCase());
						}
					}
					studyCodes.retainAll(available);
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}

				if (!studyCodes.isEmpty())
				{
					cacheStudies(studyCodes);
				}
			}
		}
		return studyCodes;
	}

	private static List<String> readStudiesFromCache()
	{
		try
		{
			if (!new File(getStudiesFileName()).exists()) return null;

			List<String> studies = new ArrayList<String>();
			BufferedReader reader = new BufferedReader(new FileReader(getStudiesFileName()));

			// Read date

			String date = reader.readLine() + "/";
			if (!broadDataURL.endsWith(date)) broadDataURL = broadDataURL.substring(0,
				broadDataURL.lastIndexOf("_") + 1) + date;

			for (String line = reader.readLine(); line != null; line = reader.readLine())
			{
				if (!line.isEmpty()) studies.add(line);
			}

			reader.close();
			return studies;
		}
		catch (IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}

	private static void cacheStudies(List<String> studies)
	{
		try
		{
			BufferedWriter writer = new BufferedWriter(
				new FileWriter(getStudiesFileName()));

			// write analysis date
			String s = getBroadDataURL();
			s = s.substring(s.indexOf("__") + 2, s.lastIndexOf("/"));
			writer.write(s + "\n");

			for (String study : studies)
			{
				writer.write(study + "\n");
			}

			writer.close();

		} catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	private static String getStudiesFileName()
	{
		return getBroadCacheDir() + CACHED_STUDIES_FILE;
	}

	private static String getBroadDateString()
	{
		String s = getBroadDataURL();
		s = s.substring(s.indexOf("__") + 2, s.lastIndexOf("/"));
		s = s.replaceAll("_", "");
		return s;
	}

	private static String getBroadDataURL(String study)
	{
		return getBroadDataURL() + "data/" + study + "/" + getBroadDateString() + "/";
	}

	private static String getBroadCacheDir()
	{
		if (cacheDir == null)
		{
			String s = BROAD_DIR;
			File f = new File(s);
			if (!f.exists()) f.mkdirs();
			return s;
		}
		return cacheDir;
	}

	private static List<String> getBroadAnalysisFileNames(String study)
	{
		List<String> list = new ArrayList<String>(30);
		try
		{
			URL url = new URL(getBroadDataURL(study));

			URLConnection con = url.openConnection();
			BufferedReader reader = new BufferedReader(new InputStreamReader(con.getInputStream()));
			for (String line = reader.readLine(); line != null; line = reader.readLine())
			{
				String start = "<li><a href=\"";
				if (line.startsWith(start))
				{
					String file = line.substring(start.length(), line.indexOf("\">"));
					list.add(file);
				}
			}
			reader.close();
		}
		catch (IOException e)
		{
			System.out.println(e);
		}
		return list;
	}

	private static String getMutsigFileName(List<String> list)
	{
		for (String s : list)
		{
			if (s.contains(MUTSIG_ANALYSIS_SUBSTR)) return s;
		}
		return null;
	}

	private static boolean mutSigAvailable(String study)
	{
		List<String> analysisFiles = getBroadAnalysisFileNames(study);
		return !analysisFiles.isEmpty() &&
			getMutsigFileName(analysisFiles) != null;
	}

	private static String getCachedMutsigFileName(String study)
	{
		return getBroadCacheDir() + study + "-mutsig.txt";
	}

	private static String getTempFileName()
	{
		return getBroadCacheDir() + "temp.tar.gz";
	}

	private static void deleteTempFile()
	{
		new File(getBroadCacheDir() + "temp.tar.gz").delete();
	}

	private static boolean downloadMutsig(String study, List<String> analysisFileNames)
	{
		String s = getMutsigFileName(analysisFileNames);
		if (s != null)
		{
			if (Download.downloadAsIs(getBroadDataURL(study) + s, getTempFileName()))
			{
				for (String name : mutsigPartialFileNames)
				{
					if (FileUtil.extractEntryContainingNameInTARGZFile(getTempFileName(), name,
						getCachedMutsigFileName(study)))
					{
						deleteTempFile();
						return true;
					}
				}
			}
		}
		return false;
	}

	public static Set<String> getMutsigGenes(String study, double thr, boolean qval)
	{
		if (!ensureStudyCached(study)) return Collections.emptySet();

		Set<String> genes = new HashSet<String>();
		genes.addAll(readGenesFromMutsig(study, thr, qval));
		return genes;
	}

	public static boolean ensureStudyCached(String study)
	{
		if (!getStudyCodes().contains(study))
		{
			System.out.println("Study " + study + " is unknown.");
			return false;
		}

		String file = getCachedMutsigFileName(study);
		if (!new File(file).exists())
		{
			downloadMutsig(study, getBroadAnalysisFileNames(study));
		}
		return new File(file).exists();
	}

	/**
	 * @param qval if true, then qval is used, else pval is used
	 */
	public static Set<String> readGenesFromMutsig(String study, double thr, boolean qval)
	{
		ensureStudyCached(study);
		Set<String> set = new HashSet<String>();
		String s = FileUtil.getFileContent(getCachedMutsigFileName(study));

		for (String line : s.split("\n"))
		{
			if (line.startsWith("rank")) continue;

			String[] token = line.split("\t");

			int index = qval ? token.length - 1 : token.length - 2;

			double val = token[index].startsWith("<") ? 0:
				Double.parseDouble(token[index]);

			if (val < thr)
			{
				String symbol = HGNC.getSymbol(token[1]);
				if (symbol != null) set.add(symbol);
			}
		}

		return set;
	}

	/**
	 * @param qval if true, then qval is used, else pval is used
	 */
	public static Map<String, Double> readGenesSignificance(String study, boolean qval)
	{
		ensureStudyCached(study);

		Map<String, Double> map = new HashMap<String, Double>();
		String s = FileUtil.getFileContent(getCachedMutsigFileName(study));

		for (String line : s.split("\n"))
		{
			if (line.startsWith("rank")) continue;

			String[] token = line.split("\t");

			int index = qval ? token.length - 1 : token.length - 2;

			double val = token[index].startsWith("<") ? 0:
				Double.parseDouble(token[index]);

			String symbol = HGNC.getSymbol(token[1]);
			if (symbol != null) map.put(symbol, val);
		}

		return map;
	}

	public static void main(String[] args) throws IOException
	{
		for (PortalDatasetEnum dataEnum : PortalDatasetEnum.values())
		{
			System.out.println("dataEnum = " + dataEnum);
			Map<String, Double> map = readGenesSignificance(dataEnum.name(), false);
			List<String> select = FDR.select(map, null, 0.1);
			System.out.println("select = " + select);
			System.out.println("\n");
		}
	}
}
