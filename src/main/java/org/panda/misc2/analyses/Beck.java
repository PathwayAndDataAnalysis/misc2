package org.panda.misc2.analyses;

import org.panda.resource.MSigDB;
import org.panda.resource.ReactomePathway;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

public class Beck
{
	public static final String DIR = "/home/ozgunbabur/Analyses/Beck/";
	public static void main(String[] args) throws IOException
	{
//		readDataAndRunEnrichment();
		zort();
	}

	public static void zort()
	{
		int a = 20;

		for (int i = 1; i < a; i++)
		{
			System.out.println("\ni = " + i);
			for (int j = 1; j < a; j++)
			{
				System.out.print("\tj = " + j + ", v = " + ((i*j) % a));
			}
		}
	}

	public static void readDataAndRunEnrichment() throws IOException
	{
		String inputFile = DIR + "table2.tsv";
		String[] header = FileUtil.readHeader(inputFile, 1);
		int regStartInd = ArrayUtil.indexOf(header, "30s ADP/Control");

		Set<String> select = new HashSet<>();
		Set<String> background = new HashSet<>();

		FileUtil.linesTabbedSkip(inputFile, 2).forEach(t ->
		{
			String gene = t[2];
			background.add(gene);
			String s = t[regStartInd] + t[regStartInd + 1] + t[regStartInd + 2];
//			if (s.contains(">") || s.contains("<"))
			if (!s.isEmpty())
			{
				select.add(gene);
			}
		});

		MSigDB go = new MSigDB();
		go.crop(name -> name.startsWith("GO_"));
		go.writeEnrichmentResults(select, background, 5, DIR + "enrichment-GO.tsv");
		MSigDB reactome = new MSigDB();
		reactome.crop(name -> name.startsWith("REACTOME_"));
		reactome.writeEnrichmentResults(select, background, 5, DIR + "enrichment-Reactome.tsv");
	}
}
