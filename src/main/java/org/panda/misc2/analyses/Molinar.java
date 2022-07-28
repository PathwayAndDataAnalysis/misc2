package org.panda.misc2.analyses;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.UniqueMaker;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

public class Molinar
{
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Molinar/";
	public static final String OUT_DIR = "/home/ozgunbabur/Analyses/Molinar/";

	public static void main(String[] args) throws IOException
	{
		convertCSVFiles();
	}

	public static void convertCSVFiles() throws IOException
	{
		for (File file : new File(DATA_DIR).listFiles())
		{
			if (!file.isDirectory() && file.getName().endsWith(".csv"))
			{
				convert(file);
			}
		}
	}

	private static void convert(File inFile) throws IOException
	{
		String outFile = OUT_DIR + inFile.getName().substring(0, inFile.getName().lastIndexOf(".")) + "/data.tsv";
		FileUtil.mkdirsOfFilePath(outFile);

		UniqueMaker uniq = new UniqueMaker();

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID\tSymbols\tFeature\tSites\tEffect\tSignedP");

		String[] header = FileUtil.readHeader(inFile.getPath());

		int symInd = ArrayUtil.indexOf(header, "Gene");
		int siteInd = ArrayUtil.indexOf(header, "Prot_Loc");
		int pInd = ArrayUtil.indexOf(header, "pvalue");
		int fcInd = ArrayUtil.indexOf(header, "log2(FC)");

		FileUtil.linesTabbedSkip1(inFile.getPath()).forEach(t ->
		{
			String sym = t[symInd];
			String sites = t[siteInd].replaceAll("\\.", "|");
			double signedP = Double.parseDouble(t[pInd]);
			if (t[fcInd].startsWith("-")) signedP = -signedP;

			String id = sym + "_" + sites.replaceAll("\\|", "_");
			id = uniq.get(id);
			FileUtil.lnwrite(id + "\t" + sym + "\tP\t" + sites + "\t\t" + signedP, writer);
		});

		writer.close();
	}
}
