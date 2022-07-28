package org.panda.misc2.analyses;

import org.panda.misc2.EdgeRReader;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;

import java.util.HashMap;
import java.util.Map;

public class Rohan
{
	public static void main(String[] args)
	{
		prepareInputs();
	}

	public static void prepareInputs()
	{
		String dataDir = "/home/ozgunbabur/Data/Rohan/";
		String outDir = "/home/ozgunbabur/Analyses/Rohan/";

		String prot = readProteomicData(dataDir + "FMD Proteomics.csv");
		String trans = readTranscriptomicData(dataDir + "FMD RNAseq-up.csv", dataDir + "FMD RNAseq-dw.csv");

		FileUtil.writeStringToFile("ID\tSymbols\tFeature\tSites\tEffect\tSignedP\n" + prot + trans, outDir + "data.tsv");
	}

	public static String readProteomicData(String filename)
	{
		String[] header = FileUtil.readHeader(filename, 5);
		int symIndex = ArrayUtil.indexOf(header, "Description");
		int fcInd = ArrayUtil.indexOf(header, "logFC");
		int pInd = ArrayUtil.indexOf(header, "PValue");

		StringBuilder builder = new StringBuilder();
		FileUtil.linesTabbed(filename).skip(6).forEach(t ->
		{
			String sym = EdgeRReader.extractGeneSymbolFromDescription(t, symIndex);
			if (sym == null) return;

			double sp = EdgeRReader.readSignedP(t, fcInd, pInd);

			builder.append(sym).append("\t").append(sym).append("\tG\t\t\t").append(sp).append("\n");
		});

		return builder.toString();
	}

	public static String readTranscriptomicData(String upFile, String dwFile)
	{
		Map<String, Double> valMap = new HashMap<>();

		FileUtil.linesTabbed(upFile).skip(3).forEach(t ->
		{
			String sym = t[1];
			if (sym.equals("NA")) return;
			if (t[3].equals("NA")) return;
			double p = Double.parseDouble(t[3]);

			valMap.put(sym, p);
		});
		FileUtil.linesTabbed(dwFile).skip(3).forEach(t ->
		{
			String sym = t[1];
			if (sym.equals("NA")) return;

			if (valMap.containsKey(sym))
			{
				System.err.println("Duplicate gene = " + sym);
			}

			if (t[3].equals("NA")) return;
			double p = Double.parseDouble(t[3]);

			valMap.put(sym, -p);
		});

		StringBuilder sb = new StringBuilder();
		for (String sym : valMap.keySet())
		{
			sb.append(sym).append("-rna\t").append(sym).append("\tR\t\t\t").append(valMap.get(sym)).append("\n");
		}
		return sb.toString();
	}
}
