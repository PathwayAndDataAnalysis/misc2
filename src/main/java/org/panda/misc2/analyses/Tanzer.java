package org.panda.misc2.analyses;

import org.panda.misc2.analyses.cptacpancan.CPTACPanCanSingleDEFileAnalysis;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * CP converter for Tanzer et al data.
 */
public class Tanzer
{
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Tanzer/";
	public static final String OUT_DIR = "/home/ozgunbabur/Analyses/Tanzer/";

	public static void main(String[] args) throws IOException
	{
		convert();
//		prepareAnalysisFolders();
	}

	public static void convert() throws IOException
	{
		String inFile = DATA_DIR + "sup-tab-1.csv";
		BufferedWriter writer = FileUtil.newBufferedWriter(OUT_DIR + "data.tsv");
		writer.write("ID\tSymbols\tSites\tFeature\tEffect");
		String[] header = FileUtil.readHeader(inFile);

		for (int i = 2; i < 9; i+=2)
		{
			String sample = header[i].substring(header[i].indexOf("q-value ") + 8).replaceAll(" ", "-");
			writer.write("\t" + sample);
		}

		FileUtil.linesTabbedSkip1(inFile).forEach(t ->
		{
			String id = t[0];
			if (id.startsWith(";")) id = id.substring(1);

			List<String> genes = Arrays.asList(id.substring(0, id.indexOf("_")).split(";"));
			genes = genes.stream().filter(g -> !g.isEmpty()).collect(Collectors.toList());
			String site = id.substring(id.indexOf("_") + 1, id.lastIndexOf("_"));

			if (site.contains("_")) return;

			if (!Character.isDigit(site.charAt(site.length() - 1)))
			{
				System.out.println("id = " + id);
				return;
			}

			if (genes.size() > 1)
			{
				String sites = site;
				for (int i = 1; i < genes.size(); i++)
				{
					sites += " " + site;
				}
				site = sites;
			}

			FileUtil.lnwrite(id + "\t" + CollectionUtil.merge(genes, " ") + "\t" + site + "\tP\t", writer);

			for (int i = 2; i < 9; i+=2)
			{
				double p = Double.parseDouble(t[i]);
				if (t[i+1].startsWith("-")) p = -p;
				FileUtil.tab_write(p, writer);
			}
		});

		writer.close();
	}

	private static void prepareAnalysisFolders() throws IOException
	{
		prepareAnalysisFoldersDiffExp(OUT_DIR + "data.tsv", OUT_DIR);
	}

	private static void prepareAnalysisFoldersDiffExp(String dataFile, String rootDir) throws IOException
	{
		FileUtil.mkdirs(rootDir);
		String dataName = dataFile.substring(dataFile.lastIndexOf("/") + 1);

		List<String> cases = List.of(FileUtil.readHeader(dataFile));
		cases = cases.subList(5, cases.size());

		for (int i = 0; i < cases.size(); i++)
		{
			String c1 = cases.get(i);
			String dir = rootDir + File.separator + c1 + File.separator;
			FileUtil.mkdirs(dir);

			String valPointers = "value-column = " + c1;

			BufferedWriter writer = new BufferedWriter(new FileWriter(dir + "parameters.txt"));
			writer.write(getParameters(dataName) + valPointers);
			writer.close();
		}
	}

	public static String getParameters(String dataName)
	{
		return "proteomics-values-file = ../data.tsv\n" +
			"id-column = ID\n" +
			"symbols-column = Symbols\n" +
			"sites-column = Sites\n" +
			"feature-column = Feature\n" +
			"effect-column = Effect\n" +
			"\n" +
			"value-transformation = signed-p-values\n" +
			"threshold-for-data-significance = 0.1 phosphoprotein\n" +
			"\n" +
			"color-saturation-value = 15\n" +
			"\n" +
			"calculate-network-significance = true\n" +
			"permutations-for-significance = 10000\n" +
			"fdr-threshold-for-network-significance = 0.1\n" +
			"use-network-significance-for-causal-reasoning = true\n" +
			"\n";
	}

}
