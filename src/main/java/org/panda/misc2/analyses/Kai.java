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

//		doGSEA(OUT_DIR + "KO-LFD-vs-WT-LFD/data.tsv", OUT_DIR + "KO-LFD-vs-WT-LFD/gsea.tsv");
		convertBlazevPaper();
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

	static void convertLinPaper() throws IOException
	{
		String base = "/home/ozgunbabur/Analyses/Kai/LinPaper/";
		StringBuilder sb = new StringBuilder();
		covertLinPhospho(base + "phosphoproteomic-data.tsv", sb);
		covertLinTotalProt(base + "proteomic-data.tsv", sb);

		String outFile = base + "/data.tsv";
		FileUtil.mkdirsOfFilePath(outFile);
		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID\tSymbols\tSites\tFeature\tEffect\t6h\t24h\t72h");
		writer.write(sb.toString());
		writer.close();
	}

	static void covertLinTotalProt(String file, StringBuilder sb)
	{
		int idInd = 2;
		int valStartInd = 11;

		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			if (t.length <= valStartInd || t[valStartInd].isEmpty()) return;

			String up = t[idInd];

			String mSym = MGI.get().getSymbol(up);
			if (mSym == null) return;

			Set<String> hSyms = MGI.get().getCorrespondingHumanSymbols(mSym);

			if (hSyms.size() == 1)
			{
				String sym = hSyms.iterator().next();

				double p6 = Double.parseDouble(t[valStartInd + 7]);
				if (t[valStartInd].startsWith("-")) p6 = -p6;

				double p24 = Double.parseDouble(t[valStartInd + 7 + 1]);
				if (t[valStartInd + 1].startsWith("-")) p24 = -p24;

				double p72 = Double.parseDouble(t[valStartInd + 7 + 2]);
				if (t[valStartInd + 2].startsWith("-")) p72 = -p72;

				sb.append("\n").append(sym).append("\t").append(sym).append("\t\tG\t\t").append(p6).append("\t").append(p24).append("\t").append(p72);
			}
		});
	}
	static void covertLinPhospho(String file, StringBuilder sb)
	{
		System.out.println("file = " + file);
		int idInd = 2;
		int valStart = 13;
		int seqInd = 4;
		int siteInd = 3;

		UniqueMaker um = new UniqueMaker();

		int[] cnt = new int[4];
		int ONE_TO_ONE = 0;
		int ONE_TO_NONE = 1;
		int ONE_TO_MANY = 2;
		int NO_HUM_SYM = 3;

		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			if (t.length < 5 || t[seqInd].isEmpty() || t[valStart].isEmpty() || t[siteInd].isEmpty()) return;

			boolean canonical = !t[idInd].contains("-");
			String up = t[idInd];
			if (!canonical) up = up.substring(0, up.indexOf("-"));

			String[] siteArr = t[siteInd].split(";");
			String[] seqArr = t[seqInd].split(";");

			List<String> siteList = new ArrayList<>();
			for (int i = 0; i < siteArr.length; i++)
			{
				String site = siteArr[i].replace("(", "").replace(")", "");

				if (!canonical)
				{
					if (seqArr[i].contains("_")) seqArr[i] = seqArr[i].substring(0, seqArr[i].indexOf("_"));
					int loc = UniProtSequence.get().getStartLocation(up, seqArr[i]);
					if (loc > 0)
					{
						loc += 15;
						siteList.add(site.substring(0, 1) + loc);
					}
				}
				else siteList.add(site);
			}

			if (siteList.isEmpty()) return; // no need to consider this row further

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

					double p6 = Double.parseDouble(t[valStart+7]);
					if (p6 == 0) p6 = 1E-10;
					if (t[valStart].startsWith("-")) p6 = -p6;

					double p24 = Double.parseDouble(t[valStart+7+1]);
					if (p24 == 0) p24 = 1E-10;
					if (t[valStart+1].startsWith("-")) p24 = -p24;

					double p72 = Double.parseDouble(t[valStart+7+2]);
					if (p72 == 0) p72 = 1E-10;
					if (t[valStart+2].startsWith("-")) p72 = -p72;

					sb.append("\n").append(id).append("\t").append(sym).append("\t")
						.append(CollectionUtil.merge(siteList, "|")).append("\tP\t\t").append(p6).append("\t").append(p24).append("\t").append(p72);
				} else cnt[NO_HUM_SYM]++;
			} else if (map.isEmpty()) cnt[ONE_TO_NONE]++;
			else cnt[ONE_TO_MANY]++;
		});
		System.out.println("cnt[ONE_TO_ONE] = " + cnt[ONE_TO_ONE]);
		System.out.println("cnt[ONE_TO_MANY] = " + cnt[ONE_TO_MANY]);
		System.out.println("cnt[ONE_TO_NONE] = " + cnt[ONE_TO_NONE]);
		System.out.println("cnt[NO_HUM_SYM] = " + cnt[NO_HUM_SYM]);
	}

	static void convertBlazevPaper() throws IOException
	{
		String base = "/home/ozgunbabur/Analyses/Kai/BlazevPaper/";
		StringBuilder sb = new StringBuilder();
		covertBlazevPhospho(base + "phosphoproteomics.tsv", sb);
		covertBlazevTotalProt(base + "proteomics.tsv", sb);

		String outFile = base + "/data.tsv";
		FileUtil.mkdirsOfFilePath(outFile);
		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID\tSymbols\tSites\tFeature\tEffect\tPostEndurance/PreEndurance\tRecoveryEndurance/PreEndurance\tPostSprint/PreSprint\tRecoverySprint/PreSprint\tPostResistance/PreResistance\tRecoveryResistance/PreResistance");
		writer.write(sb.toString());
		writer.close();
	}

	static void covertBlazevTotalProt(String file, StringBuilder sb)
	{
		int symInd = 3;
		int valStartInd = 4;

		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			if (t.length <= valStartInd || t[valStartInd].isEmpty()) return;

			if (!t[symInd].isEmpty())
			{
				String[] syms = t[symInd].split(";");

				double[] p = new double[6];

				for (int i = 0; i < p.length; i++)
				{
					p[i] = Double.parseDouble(t[valStartInd + 6 + (2 * i)]);
					if (t[valStartInd + i].startsWith("-")) p[i] = -p[i];
				}

				sb.append("\n").append(ArrayUtil.merge("-", syms)).append("\t").append(ArrayUtil.merge(" ", syms)).append("\t\tG\t");
				for (double v : p)
				{
					sb.append("\t").append(v);
				}
			}
		});
	}
	static void covertBlazevPhospho(String file, StringBuilder sb)
	{
		int idInd = 4;
		int valStart = 8;
		int siteInd = 5;

		UniqueMaker um = new UniqueMaker();

		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			if (t.length < 5 || t[valStart].isEmpty() || t[siteInd].isEmpty()) return;

			String[] ss = t[idInd].split(";");
			String sym = ss[0];
			String site = ss[1];
			if (ss.length > 2) throw new RuntimeException("More than one site or symbol?: " + t[idInd]);

			String id = sym + "-" + site;
			id = um.get(id);

			double[] p = new double[6];

			for (int i = 0; i < p.length; i++)
			{
				p[i] = Double.parseDouble(t[valStart + 6 + (3 * i)]);
				if (t[valStart + i].startsWith("-")) p[i] = -p[i];
			}

			sb.append("\n").append(id).append("\t").append(sym).append("\t").append(site).append("\tP\t");
			for (double v : p)
			{
				sb.append("\t").append(v);
			}
		});
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
