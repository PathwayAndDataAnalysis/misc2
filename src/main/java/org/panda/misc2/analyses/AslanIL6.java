package org.panda.misc2.analyses;

import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.UniqueMaker;

import java.io.BufferedWriter;
import java.io.IOException;

public class AslanIL6
{
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Aslan/IL6/";
	public static final String OUT_DIR = "/home/ozgunbabur/Analyses/Aslan-IL6/";

	public static void main(String[] args) throws IOException
	{
		convert();
	}

	public static void convert() throws IOException
	{
		String inFile = DATA_DIR + "3v3.csv";
		String[] header = FileUtil.readHeader(inFile, 5);
		int symInd = ArrayUtil.indexOf(header, "UniProt Gene Name");
//		int fdrInd = ArrayUtil.indexOf(header, "FDR");
		int fdrInd = ArrayUtil.indexOf(header, "PValue");
		int fcInd = ArrayUtil.indexOf(header, "FC");
		int siteInd = ArrayUtil.indexOf(header, "Site List");

		BufferedWriter writer = FileUtil.newBufferedWriter(OUT_DIR + "data-noadj.tsv");
		writer.write("ID\tSymbols\tFeature\tSites\tEffect\tSignedP");

		UniqueMaker um = new UniqueMaker();

		FileUtil.linesTabbedSkip(inFile, 6).forEach(t ->
		{
			String sym = t[symInd].replaceAll("\"", "");
			if (sym.contains(" ")) sym = sym.substring(0, sym.indexOf(" "));
			String sites = t[siteInd].replaceAll("; ", "|");
			double p = Double.parseDouble(t[fdrInd]);
			if (t[fcInd].startsWith("-")) p = -p;

			String id = sym + "_" + sites.replaceAll("\\|", "_");
			id = um.get(id);

			FileUtil.lnwrite(id + "\t" + sym + "\tP\t" + sites + "\t\t" + p, writer);
		});

		writer.close();
	}
}
