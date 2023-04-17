package org.panda.misc2.analyses;

import org.biopax.paxtools.model.level3.Gene;
import org.panda.misc2.KinaseEnrichment;
import org.panda.misc2.TFEnrichment;
import org.panda.resource.*;
import org.panda.utility.*;
import org.panda.utility.statistics.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class Kai
{
	public static final String[] CASES = new String[]{"KO-HFD-vs-WT-HFD", "KO-LFD-vs-WT-LFD", "WT-HFD-vs-WT-LFD"};
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Kai/v2/";
	public static final String OUT_DIR = "/home/ozgunbabur/Analyses/Kai/v2/";

	public static void main(String[] args) throws IOException
	{
//		convertAll();
//		predictActivity();

		doGSEA(OUT_DIR + "KO-LFD-vs-WT-LFD/data.tsv", OUT_DIR + "KO-LFD-vs-WT-LFD/gsea.tsv");
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
//		covertPhospho(DATA_DIR + aCase + "-phospho.csv", sb);
//		covertTotalProt(DATA_DIR + aCase + "-proteome.csv", sb);

		covertTotalProt(DATA_DIR + aCase + ".csv", sb);

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
		int ratInd = 3;
		int pInd = 4;

		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			if (t.length <= pInd || t[pInd].isEmpty() || t[ratInd].isEmpty()) return;

			String up = t[idInd];

			String mSym = MGI.get().getSymbol(up);
			if (mSym == null) return;

			Set<String> hSyms = MGI.get().getCorrespondingHumanSymbols(mSym);

			if (hSyms.size() == 1)
			{
				String sym = hSyms.iterator().next();

				double p = Double.parseDouble(t[pInd]);
				if (t[ratInd].startsWith("-")) p = -p;

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
				} else cnt[NO_HUM_SYM]++;
			} else if (map.isEmpty()) cnt[ONE_TO_NONE]++;
			else cnt[ONE_TO_MANY]++;
		});
		System.out.println("cnt[ONE_TO_ONE] = " + cnt[ONE_TO_ONE]);
		System.out.println("cnt[ONE_TO_MANY] = " + cnt[ONE_TO_MANY]);
		System.out.println("cnt[ONE_TO_NONE] = " + cnt[ONE_TO_NONE]);
		System.out.println("cnt[NO_HUM_SYM] = " + cnt[NO_HUM_SYM]);
	}

	public static void predictActivity() throws IOException
	{
		Map<String, Map<String, Map<String, Boolean>>> rawKSPriors = KinaseEnrichment.readPriors();
		Map<String, Map<String, Boolean>> rawTFTPriors = TFEnrichment.readPriors();

		for (File dir : new File(OUT_DIR).listFiles())
		{
			System.out.println("dir = " + dir);
			String dataFile = dir.getPath() + "/data.tsv";
			List<String> list = KinaseEnrichment.readRankedIDsFromCPFile(dataFile, "SignedP");
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

	public static void doGSEA(String inputCPFile, String outFile) throws IOException
	{
		String[] header = FileUtil.readHeader(inputCPFile);
		int symInd = ArrayUtil.indexOf(header, "Symbols");
		int pInd = ArrayUtil.indexOf(header, "SignedP");
		Map<String, Double> pMap = FileUtil.linesTabbedSkip1(inputCPFile)
			.collect(Collectors.toMap(t -> t[symInd], t -> Math.abs(Double.parseDouble(t[pInd])), Math::min));

		List<String> select = FDR.select(pMap, null, 0.1);
		Set<String> background = pMap.keySet();

//		Map<String, Set<String>> reactomeWithID = ReactomePathway.get().getAllPathways();
//		Map<String, Set<String>> reactomeWithName = reactomeWithID.keySet().stream().collect(Collectors.toMap(id -> ReactomePathway.get().getName(id), reactomeWithID::get, (strings, strings2) -> strings));
//		EnrichmentAnalyzer.run(reactomeWithName, new HashSet<>(select), background, 3, 0.1, outFile);

//		Map<String, Set<String>> reactome = MSigDB.get().getSetsNameFiltered(name -> name.startsWith("REACTOME"));
		Map<String, Set<String>> reactome = MSigDB.get().getSetsNameFiltered(name -> true);
		EnrichmentAnalyzer.run(reactome, new HashSet<>(select), background, 3, 0.1, outFile);
	}



}
