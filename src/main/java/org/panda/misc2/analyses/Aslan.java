package org.panda.misc2.analyses;

import org.panda.resource.ChEBI;
import org.panda.resource.GWASCatalog;
import org.panda.utility.ArrayUtil;
import org.panda.utility.EnrichmentAnalyzer;
import org.panda.utility.FileUtil;
import org.panda.utility.SIFFileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class Aslan
{
	public static void main(String[] args) throws IOException
	{
//		doGWASEnrichment();
		metabolicNetworkLookup();
	}

	public static void doGWASEnrichment() throws IOException
	{
		String dir = "/home/ozgunbabur/Analyses/Platelet-Blood-paper/cond2/";
		Set<String>[] selAndBg = readSignifcantGenesAndBackground(dir + "data-fdr0.1.txt",
			"Symbols", "Fold change", 0.1);
		Map<String, Set<String>> geneSets = GWASCatalog.get().getGeneSets();
		EnrichmentAnalyzer.run(geneSets, selAndBg[0], selAndBg[1], 5, 1, dir + "gwas-enrichment.tsv");
	}

	private static Set<String>[] readSignifcantGenesAndBackground(String cpDataFile, String geneColumn, String valColumn, double thr)
	{
		Set<String> background = new HashSet<>();
		Set<String> selected = new HashSet<>();
		String[] header = FileUtil.readHeader(cpDataFile);
		int geneInd = ArrayUtil.indexOf(header, geneColumn);
		int valInd = ArrayUtil.indexOf(header, valColumn);

		FileUtil.linesTabbedSkip1(cpDataFile).forEach(t ->
		{
			String gene = t[geneInd];
			background.add(gene);

			double val = Double.parseDouble(t[valInd]);
			if (Math.abs(val) <= thr)
			{
				selected.add(gene);
			}
		});

		return new Set[]{selected, background};
	}

	private static Set<String> readChEBIIDs()
	{
		return FileUtil.linesTabbedSkip1("/home/ozgunbabur/Data/Aslan/oxylipindata.csv")
			.filter(t -> !t[2].equals("N/A")).map(t -> "CHEBI:" + t[2]).collect(Collectors.toSet());
	}

	private static void metabolicNetworkLookup() throws IOException
	{
		String outDir = "/home/ozgunbabur/Data/Aslan/";
		Set<String> chebi = readChEBIIDs();
		String priorFile = "/home/ozgunbabur/Data/causal-priors.txt";
		String outFile = outDir + "metabolic-paths-between.sif";
		SIFFileUtil.writePatsBetween(priorFile, chebi, outFile);
		replaceChEBIIDsWithNames(outFile);

		String nDir = outDir + "metabolic-neighborhoods/";
		FileUtil.mkdirs(nDir);

		for (String id : chebi)
		{
			String name = ChEBI.get().getName(id);
			outFile = nDir + name.replaceAll(" ", "-") + ".sif";
			SIFFileUtil.writeNeighborhood(priorFile, Collections.singleton(id), outFile);
			replaceChEBIIDsWithNames(outFile);
		}
	}

	public static void replaceChEBIIDsWithNames(String sifFile) throws IOException
	{
		Set<String> rels = FileUtil.lines(sifFile).collect(Collectors.toSet());
		String prefix = "CHEBI:";
		String noname = "noname";
		BufferedWriter writer = FileUtil.newBufferedWriter(sifFile);
		for (String rel : rels)
		{
			String[] t = rel.split("\t");
			if (t.length < 3) continue;

			if (t[0].startsWith(prefix)) t[0] = ChEBI.get().getName(t[0]);
			if (t[0] == null) t[0] = noname;
			if (t[2].startsWith(prefix)) t[2] = ChEBI.get().getName(t[2]);
			if (t[2] == null) t[2] = noname;
			FileUtil.writeln(ArrayUtil.merge("\t", t), writer);
		}
		writer.close();
	}

}
