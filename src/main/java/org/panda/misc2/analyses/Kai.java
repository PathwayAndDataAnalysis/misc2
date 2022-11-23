package org.panda.misc2.analyses;

import org.panda.misc2.KinaseEnrichment;
import org.panda.misc2.TFEnrichment;
import org.panda.resource.HGNC;
import org.panda.resource.MGI;
import org.panda.resource.SiteMappingMouseToHuman;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.UniqueMaker;
import org.panda.utility.statistics.RankedListSignedEnrichment;
import org.panda.utility.statistics.RankedListSignedGroupedDifferentialEnrichment;
import org.panda.utility.statistics.RankedListSignedGroupedEnrichment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class Kai
{
	public static final String[] CASES = new String[]{"LFD-ins-vs-noins", "HFD-ins-vs-noins", "HFD-vs-LFD-ins", "HFD-vs-LFD-noins"};
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Kai/";
	public static final String OUT_DIR = "/home/ozgunbabur/Analyses/Kai/";

	public static void main(String[] args) throws IOException
	{
//		convertAll();
		predictActivity();
	}

	static void convertAll() throws IOException
	{
		for (String aCase : CASES)
		{
			convertCase(aCase);
		}
	}

	static void convertCase(String aCase) throws IOException
	{
		StringBuilder sb = new StringBuilder();
		covertPhospho(DATA_DIR + aCase + "-phospho.csv", sb);
		covertTotalProt(DATA_DIR + aCase + "-proteome.csv", sb);

		String outFile = OUT_DIR + aCase + "/data.tsv";
		FileUtil.mkdirsOfFilePath(outFile);
		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID\tSymbols\tSites\tFeature\tEffect\tSignedP");
		writer.write(sb.toString());
		writer.close();
	}

	static void covertTotalProt(String file, StringBuilder sb)
	{
		int idInd = 0;
		int ratInd = 1;
		int pInd = 2;

		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			if (t.length < 3 || t[pInd].isEmpty() || t[ratInd].isEmpty()) return;

			String up = t[idInd];

			String mSym = MGI.get().getSymbol(up);
			if (mSym == null) return;

			Set<String> hSyms = MGI.get().getCorrespondingHumanSymbols(mSym);

			if (hSyms.size() == 1)
			{
				String sym = hSyms.iterator().next();

				double p = Double.parseDouble(t[pInd]);
				if (t[ratInd].startsWith("0.")) p = -p;

				sb.append("\n").append(sym).append("\t").append(sym).append("\t\tG\t\t").append(p);
			}
		});
	}

	static void covertPhospho(String file, StringBuilder sb)
	{
		System.out.println("file = " + file);
		int idInd = 0;
		int ratInd = 1;
		int pInd = 2;
		int siteInd = 4;

		UniqueMaker um = new UniqueMaker();

		int[] cnt = new int[4];
		int ONE_TO_ONE = 0;
		int ONE_TO_NONE = 1;
		int ONE_TO_MANY = 2;
		int NO_HUM_SYM = 3;

		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			if (t.length < 5 || t[pInd].isEmpty() || t[ratInd].isEmpty() || t[siteInd].isEmpty()) return;

			String up = t[idInd];

			String s = t[siteInd].substring(t[siteInd].indexOf("[") + 1, t[siteInd].lastIndexOf("]"));

			List<String> siteList = new ArrayList<>();
			for (String ss : s.split("; "))
			{
				String site = ss.substring(0, ss.indexOf("("));
				siteList.add(site);
			}

			Map<String, List<String>> map = SiteMappingMouseToHuman.get().mapToHumanSite(up, siteList.toArray(new String[0]));


			if (map.size() == 1)
			{
				String hUP = map.keySet().iterator().next();
				String sym = HGNC.get().getSymbol(hUP);

				if (sym != null)
				{
					cnt[ONE_TO_ONE]++;
					siteList = map.get(hUP);

					String id = sym + "-" + CollectionUtil.merge(siteList, "-");
					id = um.get(id);

					double p = Double.parseDouble(t[pInd]);
					if (t[ratInd].startsWith("0.")) p = -p;

					sb.append("\n").append(id).append("\t").append(sym).append("\t")
						.append(CollectionUtil.merge(siteList, "|")).append("\tP\t\t").append(p);
				}
				else cnt[NO_HUM_SYM]++;
			}
			else if (map.isEmpty()) cnt[ONE_TO_NONE]++;
			else cnt[ONE_TO_MANY]++;
		});
		System.out.println("cnt[ONE_TO_ONE] = " + cnt[ONE_TO_ONE]);
		System.out.println("cnt[ONE_TO_MANY] = " + cnt[ONE_TO_MANY]);
		System.out.println("cnt[ONE_TO_NONE] = " + cnt[ONE_TO_NONE]);
		System.out.println("cnt[NO_HUM_SYM] = " + cnt[NO_HUM_SYM]);
	}

	public static void predictActivity() throws IOException {
		Map<String, Map<String, Map<String, Boolean>>> rawKSPriors = KinaseEnrichment.readPriors();
		Map<String, Map<String, Boolean>> rawTFTPriors = TFEnrichment.readPriors();

		for (File dir : new File(OUT_DIR).listFiles())
		{
			System.out.println("dir = " + dir);
			String dataFile = dir.getPath() + "/data.tsv";
			List<String> list = KinaseEnrichment.readRankedIDsFromCPFile(dataFile);
			Map<String, Map<String, Set<String>>> geneToSiteToID = KinaseEnrichment.readDataMappingFromCPFile(dataFile);
			Map<String, Set<Map<String, Boolean>>> priors = KinaseEnrichment.convertPriors(geneToSiteToID, rawKSPriors, list);
			if (!priors.isEmpty())
			{
				RankedListSignedGroupedEnrichment.reportEnrichment(list, priors, 1000000,
					dir.getPath() + "/kinase-enrichment.tsv");
			}

			list = TFEnrichment.readRankedIDsFromCPFile(dataFile, "G");
			Map<String, Map<String, Boolean>> priors2 = TFEnrichment.convertPriors(rawTFTPriors, list, 3);
			if (!priors2.isEmpty())
			{
				RankedListSignedEnrichment.reportEnrichment(list, priors2, 1000000,
					dir.getPath() + "/tf-enrichment.tsv");
			}
		}
	}

}
