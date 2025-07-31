package org.panda.misc2.analyses;

import org.panda.resource.ChEBI;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.statistics.Summary;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class CPTACLSCC
{
	public static final String DATA_BASE = "/home/ozgunbabur/Data/CPTAC-LSCC/v3/";
	public static final String OUT_BASE = DATA_BASE.replace("Data", "Analyses");

	public static void main(String[] args) throws IOException
	{
//		prepareData();
//		prepareAnalysisDirsSubtypeAgainstOthers();
//		prepareAnalysisDirsTumorsVsNormals();
//		extractMetabolicNetwork();
//		printResultOverlap(OUT_BASE + "tumors-vs-normals-per-NMF-subtype");

//		generateNeighborhoodGraphs(OUT_BASE + "altered-subset-vs-normals", "NRF2-targets", NRF2_TARGETS);
//		generateNeighborhoodGraphs(OUT_BASE + "tumors-vs-normals-per-NMF-subtype", "keap1-nrf2-neigh",
//			new HashSet<>(Arrays.asList("KEAP1", "NFE2L2")));

//		generateWithinGSGraphs(OUT_BASE + "tumors-vs-normals-per-NMF-subtype", "keap1-nrf2-within-geneset", KEAP1GeneSet);

//		prepareDirForAlteredSubsetVsNormals(Arrays.asList("KEAP1", "NFE2L2", "CUL3"), Arrays.asList("-1", "1", "-1"));

//		printNRF2Overlap();

		Map<String, String> map = getMetabIDToNameMap();
		String siffile = OUT_BASE + "altered-subset-vs-normals/KEAP1-NFE2L2-CUL3-altered/NFE2L2-downstream-2.sif";
		String fmtfile = OUT_BASE + "altered-subset-vs-normals/KEAP1-NFE2L2-CUL3-altered/NFE2L2-downstream-2.format";
		replaceCHEBIWithNames(map, siffile, siffile, true);
		replaceCHEBIWithNames(map, fmtfile, fmtfile, true);
	}

	public static void prepareAnalysisDirsSubtypeAgainstOthers()
	{
		Map<String, List<String>> sMap = loadSubtypeMap();
		List<String> samples = loadSamples();
		String base = OUT_BASE + "NMF-subtype-vs-others";
		for (String subtype : sMap.keySet())
		{
			String dir = base + "/" + subtype + "/";
			FileUtil.mkdirs(dir);

			BufferedWriter writer = FileUtil.newBufferedWriter(dir + "parameters.txt");
			FileUtil.write(PARAM_PREFIX, writer);

			List<String> includeList = sMap.get(subtype);
			for (String sample : includeList)
			{
				FileUtil.lnwrite("test-value-column = " + sample, writer);
			}
			for (String sample : samples)
			{
				if (!sample.endsWith(".N") && !includeList.contains(sample))
				{
					FileUtil.lnwrite("control-value-column = " + sample, writer);
				}
			}

			FileUtil.closeWriter(writer);
		}
	}

	public static void prepareAnalysisDirsTumorsVsNormals()
	{
		Map<String, List<String>> sMap = loadSubtypeMap();
		List<String> samples = loadSamples();
		String base = OUT_BASE + "tumors-vs-normals-per-NMF-subtype";

		Set<String> cancers = new HashSet<>();
		Set<String> normals = new HashSet<>();

		for (String subtype : sMap.keySet())
		{
			String dir = base + "/" + subtype + "/";
			FileUtil.mkdirs(dir);

			BufferedWriter writer = FileUtil.newBufferedWriter(dir + "parameters.txt");
			FileUtil.write(PARAM_PREFIX, writer);

			List<String> includeList = sMap.get(subtype);

			for (String sample : includeList)
			{
				cancers.add(sample);
				FileUtil.lnwrite("test-value-column = " + sample, writer);
			}
			for (String sample : samples)
			{
				if (sample.endsWith(".N") && includeList.contains(sample.substring(0, sample.length() - 2)))
				{
					normals.add(sample);
					FileUtil.lnwrite("control-value-column = " + sample, writer);
				}
			}

			FileUtil.closeWriter(writer);
		}

		String dir = base + "/without-subtypes/";
		FileUtil.mkdirs(dir);

		BufferedWriter writer = FileUtil.newBufferedWriter(dir + "parameters.txt");
		FileUtil.write(PARAM_PREFIX, writer);

		for (String sample : cancers)
		{
			FileUtil.lnwrite("test-value-column = " + sample, writer);
		}
		for (String sample : normals)
		{
			FileUtil.lnwrite("control-value-column = " + sample, writer);
		}
		FileUtil.closeWriter(writer);
	}

	private static List<String> loadSamples()
	{
		return FileUtil.lines(DATA_BASE + "LSCC-combined-v3-sample-annotation.csv").skip(1)
			.map(l -> l.substring(0, l.indexOf(","))).collect(Collectors.toList());
	}

	private static Map<String, List<String>> loadSubtypeMap()
	{
		String[] header = FileUtil.lines(DATA_BASE + "LSCC-combined-v3-sample-annotation.csv").findFirst().get().split(",");
		int subtypeInd = ArrayUtil.indexOf(header, "NMF.consensus");

		Map<String, List<String>> map = new HashMap<>();

		FileUtil.lines(DATA_BASE + "LSCC-combined-v3-sample-annotation.csv").skip(1)
			.map(l -> l.split(",")).filter(t -> t.length > subtypeInd && !t[subtypeInd].isEmpty()).forEach(t ->
			{
				if (!map.containsKey(t[subtypeInd])) map.put(t[subtypeInd], new ArrayList<>());
				map.get(t[subtypeInd]).add(t[0]);
			});

		return map;
	}

	public static void prepareData()
	{
		List<String> samples = loadSamples();

		BufferedWriter writer = FileUtil.newBufferedWriter(OUT_BASE + "data.tsv");
		FileUtil.write("ID\tSymbols\tSites\tEffect\tFeature", writer);
		samples.stream().forEach(s -> FileUtil.tab_write(s, writer));

		String protFile = DATA_BASE + "LSCC-combined-v3-proteome-SpectrumMill-ratio-bridging-batch-correct-QCfilter-NArm.gct";
		String phosphoFile = DATA_BASE + "LSCC-combined-v1-phosphoproteome-SpectrumMill-ratio-bridging-batch-correct-QCfilter-NArm.gct";
		String acetylFile = DATA_BASE + "LSCC-combined-v3-acetylome-SpectrumMill-ratio-bridging-batch-correct-QCfilter-NArm.gct";
		String rnaFile = DATA_BASE + "LSCC-combined-v3-rnaseq-log2-TPM-median-norm-QCfilter-NArm.gct";
		String metabolFile = DATA_BASE + "LSCC-combined-v3-metabolomics-log2-median-norm.gct";

		processGCT(protFile, false, false, null,"G", samples, writer);
		processGCT(phosphoFile, true,false, null, "P", samples, writer);
		processGCT(acetylFile, true,false,null,"A", samples, writer);
		processGCT(rnaFile, false, false, null,"R", samples, writer);
		processGCT(metabolFile, false, true, loadChEBIMap(), "C", samples, writer);

		FileUtil.closeWriter(writer);
	}

	private static Map<String, String> loadChEBIMap()
	{
		return FileUtil.linesTabbed(DATA_BASE + "metabolite_ID_map2.csv")
			.filter(t -> !t[0].equals("NA") && !t[2].equals("NA"))
			.collect(Collectors.toMap(t -> t[0], t -> "CHEBI:" + t[2]));
	}


	private static void processGCT(String filename, boolean siteSpec, boolean metabolic, Map<String, String> chebiMap, String feature, List<String> samples, BufferedWriter writer)
	{
		String[] stats = FileUtil.lines(filename).skip(1).findFirst().get().split("\t");
		String[] header = FileUtil.lines(filename).skip(2).findFirst().get().split("\t");

		Map<String, Integer> sampleLocMap = samples.stream().collect(Collectors.toMap(s -> s, s -> ArrayUtil.indexOf(header, s)));

		int skip = Integer.parseInt(stats[3]) + 3;

		int symInd = ArrayUtil.indexOf(header, metabolic ? "HMDB_ID" : "geneSymbol");

		FileUtil.linesTabbedSkip(filename, skip).filter(t -> !t[symInd].isEmpty()).forEach(t ->
		{
			String sym = t[symInd];

			if (metabolic)
			{
				if (!chebiMap.containsKey(sym)) return;

				sym = chebiMap.get(sym);
			}

			String id;
			String sitesStr = "";

			if (siteSpec)
			{
				List<String> sites = parseSites(t[0]);
				id = sym;
				for (String site : sites)
				{
					id += "_" + site;
					sitesStr += "|" + site;
				}
				sitesStr = sitesStr.substring(1);
			}
			else
			{
				id = t[0];
			}
			id += "_" + feature;

			FileUtil.lnwrite(id + "\t" + sym + "\t" + sitesStr + "\t\t" + feature, writer);

			for (String sample : samples)
			{
				int index = sampleLocMap.get(sample);
				String val = index < 0 || t[index].equals("NA") ? "NaN" : t[index];
				FileUtil.tab_write(val, writer);
			}
		});
	}

	private static List<String> parseSites(String idInFile)
	{
		List<String> list = new ArrayList<>();
		for (String s : idInFile.split("_")[1].split("[a-z]"))
		{
			if (!s.isEmpty())
			{
				String aa = s.substring(0, 1);
				int site = Integer.parseInt(s.substring(1));

				list.add(s);
			}
		}
		return list;
	}

	public static final String PARAM_PREFIX = "proteomics-values-file = ../../data.tsv\n" +
		"id-column = ID\n" +
		"symbols-column = Symbols\n" +
		"sites-column = Sites\n" +
		"feature-column = Feature\n" +
		"effect-column = Effect\n" +
		"\n" +
		"value-transformation = significant-change-of-mean\n" +
		"fdr-threshold-for-data-significance = 0.1 protein\n" +
		"fdr-threshold-for-data-significance = 0.1 phosphoprotein\n" +
		"fdr-threshold-for-data-significance = 0.1 acetylprotein\n" +
		"fdr-threshold-for-data-significance = 0.1 rna\n" +
		"fdr-threshold-for-data-significance = 0.1 metabolite\n" +
		"\n" +
		"color-saturation-value = 15\n" +
		"\n" +
		"calculate-network-significance = true\n" +
		"permutations-for-significance = 10000\n" +
		"fdr-threshold-for-network-significance = 0.1\n" +
		"use-network-significance-for-causal-reasoning = true\n" +
		"\n" +
		"show-all-genes-with-proteomic-data = true\n" +
		"\n" +
		"data-type-for-expressional-targets = rna\n";

	private static Map<String, String> getMetabIDToNameMap()
	{
		Map<String, String> chEBIMap = loadChEBIMap();
		Map<String, String> hmToName = FileUtil.linesTabbed(DATA_BASE + "LSCC-combined-v3-metabolomics-log2-median-norm.gct")
			.filter(t -> t.length > 7 && t[5].startsWith("HMDB") && !t[7].equals("NA")).collect(Collectors.toMap(t -> t[5], t -> t[7], (v1, v2) -> v1));

		Map<String, String> map = new HashMap<>();
		hmToName.forEach((hm, name) ->
		{
			if (chEBIMap.containsKey(hm)) map.put(chEBIMap.get(hm), name);
		});
		return map;
	}

	private static void extractMetabolicNetwork()
	{
		Map<String, String> mToName = getMetabIDToNameMap();

//		String root = "NMF-subtype-vs-others";
		String root = "tumors-vs-normals-per-NMF-subtype";
		for (File dir : new File(OUT_BASE + root).listFiles())
		{
			if (FileUtil.exists(dir + "/causative.sif"))
			{
				replaceCHEBIWithNames(mToName, dir + "/causative.sif", dir + "/metabolic.sif", false);
				replaceCHEBIWithNames(mToName, dir + "/causative.format", dir + "/metabolic.format", true);
			}
		}
	}

	private static void replaceCHEBIWithNames(Map<String, String> mToName, String inputFile, String outFile, boolean keepOtherLines)
	{
		List<String> lines = new ArrayList<>();
		FileUtil.lines(inputFile).filter(l -> l.contains("\t")).forEach(l ->
		{
			boolean hasChEBI = false;
			while (l.contains("CHEBI:"))
			{
				hasChEBI = true;
				int ind = l.indexOf("CHEBI:");
				int end = l.indexOf("\t", ind + 6);

				if (end < 0) System.out.println(l);

				String chID = l.substring(ind, end);

				String name = ChEBI.get().getName(chID);
				if (name == null) name = mToName.get(chID);
				if (name == null) name = "noname";

				l = l.replace(chID, name);
			}
			if (hasChEBI || keepOtherLines)
			{
				lines.add(l);
			}
		});

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		lines.forEach(l -> FileUtil.writeln(l, writer));
		FileUtil.closeWriter(writer);
	}

	private static void printResultOverlap(String dir)
	{
		Map<String, Set<String>> relMap = new HashMap<>();

		for (File subdir : new File(dir).listFiles())
		{
			Set<String> rels = FileUtil.linesTabbed(subdir.getPath() + "/causative.sif")
				.filter(t -> t.length > 2)
				.map(t -> t[0] + "\t" + t[1] + "\t" + t[2])
				.collect(Collectors.toSet());

			relMap.put(subdir.getName(), rels);
		}

		String[] typeNames = new String[relMap.size()];
		Collection[] cols = new Collection[relMap.size()];

		int i = 0;
		for (String name : relMap.keySet())
		{
			typeNames[i] = name;
			cols[i++] = relMap.get(name);
		}

		CollectionUtil.printNameMapping(typeNames);
		CollectionUtil.printVennCounts(cols);
	}



	private static void generateNeighborhoodGraphs(String baseDir, String outName, Set<String> goi) throws IOException
	{
		for (File dir : new File(baseDir).listFiles())
		{
			if (FileUtil.exists(dir + "/causative.sif"))
			{
				generateNeighborhoodGraph(dir + "/causative.sif", dir + "/" + outName + ".sif", goi);
			}
		}
	}

	private static void generateWithinGSGraphs(String baseDir, String outName, Set<String> goi) throws IOException
	{
		for (File dir : new File(baseDir).listFiles())
		{
			if (FileUtil.exists(dir + "/causative.sif"))
			{
				generateWithinGSGraph(dir + "/causative.sif", dir + "/" + outName + ".sif", goi);
			}
		}
	}

	private static void generateNeighborhoodGraph(String inSIFFile, String outSIFFile, Set<String> goi) throws IOException
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(outSIFFile);

		FileUtil.lines(inSIFFile).forEach(l ->
		{
			String[] t = l.split("\t");
			if (t.length > 2)
			{
				if (goi.contains(t[0]) || goi.contains(t[2])) FileUtil.writeln(l, writer);
			}
		});

		FileUtil.closeWriter(writer);
		String inFmt = inSIFFile.substring(0, inSIFFile.lastIndexOf(".")) + ".format";
		String outFmt = outSIFFile.substring(0, outSIFFile.lastIndexOf(".")) + ".format";
		FileUtil.copyFile(inFmt, outFmt);

		Map<String, String> mToName = getMetabIDToNameMap();
		replaceCHEBIWithNames(mToName, outSIFFile, outSIFFile, true);
		replaceCHEBIWithNames(mToName, outFmt, outFmt, true);
	}

	private static void generateWithinGSGraph(String inSIFFile, String outSIFFile, Set<String> goi) throws IOException
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(outSIFFile);

		FileUtil.lines(inSIFFile).forEach(l ->
		{
			String[] t = l.split("\t");
			if (t.length > 2)
			{
				if (goi.contains(t[0]) && goi.contains(t[2])) FileUtil.writeln(l, writer);
			}
		});

		FileUtil.closeWriter(writer);
		String inFmt = inSIFFile.substring(0, inSIFFile.lastIndexOf(".")) + ".format";
		String outFmt = outSIFFile.substring(0, outSIFFile.lastIndexOf(".")) + ".format";
		FileUtil.copyFile(inFmt, outFmt);
	}

	private static final Set<String> KEAP1GeneSet = new HashSet<>(Arrays.asList("ME1 GCLC MAFG ABCC1 NQO1 TXN NFE2L2 PSMD11 SLC7A11 KEAP1 TALDO1 SRXN1 TXNRD1 PGD BRCA1 TKT AMER1 CUL1 MYC GSR GCLM G6PD MUL1 GSTA1 PRKCI ABCC3 NFKB1 EGF UBXN7 NPLOC4 CCL2 PRDX1 ABCF2 PSMB7 PRKAA2 IDH1 SEM1 SQSTM1 SESN1 SESN2 CSNK2A2 PSMD3 PSMD2 VCP PSMB6 UBA52 PSMD12 PSMB2 PSMB4 RPS27A UFD1 PALB2 PSMB5 MAP1LC3B PSMC5 PSMD14 SP1 CSNK2A1 ATF4 PSMD6 AREG PSMD8 BCL2L1 UBB PSMD7 NOTCH1 PSMB3 PSMC1 PDGFA TRIM21 PSMC2 SOD3 PSMA5 CDKN1A CREBBP CDKN2A RBX1 PSMA7 PSMD13 FBXL17 PSMA2 BACH1 PSMA4 PSMC4 ADRM1 AKT3 PRKCD BCL2 EP300 RELA GSK3B CSNK2B SKP2 CUL3 HMOX1 PSMA3 BTRC PSMA1 PSMC3 EIF2AK3 PSMD1 ABCG2 AKT2 PSMB1 DPP3 UBC SKP1 AKT1 GSTA3 PSMC6 PSMA6 MAFK".split(" ")));

	private static List<String> getSamplesAltered(List<String> genes, List<String> expectedCNV)
	{
		String[] header = FileUtil.lines(DATA_BASE + "LSCC-combined-v3-sample-annotation.csv").findFirst().get().split(",");

		List<Integer> mutInds = new ArrayList<>(genes.size());
		List<Integer> cnvInds = new ArrayList<>(genes.size());
		List<String> cnvVals = new ArrayList<>(genes.size());

		for (int i = 0; i < genes.size(); i++)
		{
			String gene = genes.get(i);
			int index = ArrayUtil.indexOf(header, gene + ".mutation.status");
			if (index > 0) mutInds.add(index);
			index = ArrayUtil.indexOf(header, gene + ".CNV");
			if (index > 0)
			{
				cnvInds.add(index);
				cnvVals.add(expectedCNV.get(i));
			}
		}

		int maxCol = Summary.max(cnvInds);
		List<String> samples = new ArrayList<>();

		FileUtil.lines(DATA_BASE + "LSCC-combined-v3-sample-annotation.csv").skip(1)
			.map(l -> l.split(",")).filter(t -> t.length > maxCol && !t[0].endsWith(".N")).forEach(t ->
			{
				boolean hit = false;
				for (int mutInd : mutInds)
				{
					if (t[mutInd].equals("1"))
					{
						hit = true;
						break;
					}
				}
				if (!hit)
				{
					for (int i = 0; i < cnvInds.size(); i++)
					{
						int cnvInd = cnvInds.get(i);
						if (t[cnvInd].equals(cnvVals.get(i)))
						{
							hit = true;
							break;
						}
					}
				}
				if (hit) samples.add(t[0]);
			});

		return samples;
	}

	private static List<String> getCorrespondingNormals(List<String> soi)
	{
		List<String> samples = loadSamples();
		List<String> normals = new ArrayList<>(soi.size());
		for (String sample : samples)
		{
			if (sample.endsWith(".N"))
			{
				String id = sample.substring(0, sample.length()-2);
				if (soi.contains(id)) normals.add(sample);
			}
		}
		return normals;
	}

	private static void prepareDirForAlteredSubsetVsNormals(List<String> genes, List<String> expectedCNV)
	{
		String outBase = OUT_BASE + "altered-subset-vs-normals/";
		String dir = "";
		for (String gene : genes)
		{
			dir += gene + "-";
		}
		dir += "altered/";
		dir = outBase + dir;

		FileUtil.mkdirs(dir);

		BufferedWriter writer = FileUtil.newBufferedWriter(dir + "parameters.txt");
		FileUtil.writeln(PARAM_PREFIX, writer);
		List<String> samplesAltered = getSamplesAltered(genes, expectedCNV);
		samplesAltered.forEach(s -> FileUtil.lnwrite("test-value-column = " + s, writer));
		getCorrespondingNormals(samplesAltered).forEach(s -> FileUtil.lnwrite("control-value-column = " + s, writer));
		FileUtil.closeWriter(writer);
	}

	public static final Set<String> NRF2_TARGETS = new HashSet<>(Arrays.asList("G6PD GCLM GPX2 SRXN1 TXN UGT1A7 ALDH3A1 CBR1 GCLC GSR ME1 NQO1 PGD PTGR1 SLC7A11 TXNRD1CBR3 CES1 CYP4F3 EPHX1 TSPAN7 ABCC1 AKR1B10 AKR1C1 AKR1C3 CYP4F11 GSTM3 MAP2 SPP1 UGT1A6 GSTA1 ADAM23 ADH7 AKR1C2 SOST WNT5AEPS8 MAP1B ABCB6 AKR1C4 NR0B1 UCHL1 WNT11 NRCAM CBX2 COCH GSTM2 RAB3B RAB6B SCIN SOX2 ALDH1A1 CARD11 GPAT3 KRT7 NGFR S100A9CA9 CAPN13 COL21A1 CXCL6 CYP4X1 EDARADD FBXO2 MMP13 PADI3 SULT1E1 TCN1".split(" ")));

	private static void printNRF2Overlap()
	{
		Set<String> cpTargs = FileUtil.linesTabbed("/home/ozgunbabur/Documents/causal-priors.txt")
			.filter(t -> t.length > 2 && t[0].equals("NFE2L2") && (t[1].startsWith("upr") || t[1].startsWith("downr")))
			.map(t -> t[2]).collect(Collectors.toSet());

		CollectionUtil.printVennSets(cpTargs, NRF2_TARGETS);

		System.out.println("\n\n");

		CollectionUtil.printVennSets(KEAP1GeneSet, NRF2_TARGETS);
	}
}
