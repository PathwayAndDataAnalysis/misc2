package org.panda.misc2.analyses;

import org.panda.resource.HGNC;
import org.panda.resource.MGI;
import org.panda.resource.SiteMappingMouseToHuman;
import org.panda.resource.network.UniProt;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.UniqueMaker;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Newton
{
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Newton/";
	public static final String OUT_DIR = "/home/ozgunbabur/Analyses/Newton/";

	final static String[] valColNames = new String[]{"WT_Rep1","WT_Rep2", "WT_Rep3","HET_Rep1","HET_Rep2","HET_Rep3","HOM_Rep1","HOM_Rep2","HOM_Rep3"};

	public static void main(String[] args) throws IOException
	{
		convert();
	}

	public static void convert() throws IOException
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(OUT_DIR + "data.tsv");
		writer.write("ID\tSymbols\tSites\tFeature\tEffect\tSignedP");
//		for (String valColName : valColNames)
//		{
//			writer.write("\t" + valColName);
//		}

		writer.write(convertPhospho(DATA_DIR + "phosphoproteomics.csv"));
//		writer.write(convertTotProt(DATA_DIR + "globalproteomics.csv"));
		writer.close();
	}


	public static String convertPhospho(String inFile) {
		System.out.println("Newton.convertPhospho");
		StringBuilder sb = new StringBuilder();

		String[] header = FileUtil.readHeader(inFile);
		int upInd = ArrayUtil.indexOf(header, "Protein");
		int siteInd = ArrayUtil.indexOf(header, "ProtLoc");
		int pInd = ArrayUtil.indexOf(header, "WTvHOM_pvalue");
		int fcInd = ArrayUtil.indexOf(header, "WTvHOM_FC");


		UniqueMaker um = new UniqueMaker();

		int[] cnt = new int[4];
		int ONE_TO_ONE = 0;
		int ONE_TO_NONE = 1;
		int ONE_TO_MANY = 2;
		int NO_HUM_SYM = 3;

		FileUtil.linesTabbedSkip1(inFile).forEach(t ->
		{
			String mup = t[upInd];

			String[] site = t[siteInd].split("\\.");

			Map<String, List<String>> mappedSites = SiteMappingMouseToHuman.get().mapToHumanSite(mup, site);

			if (mappedSites.size() == 1)
			{
				cnt[ONE_TO_ONE]++;

				String up = mappedSites.keySet().iterator().next();
				List<String> sites = mappedSites.get(up);

				String sym = HGNC.get().getSymbol(up);
				if (sym == null)
				{
					cnt[NO_HUM_SYM]++;
					return;
				}

				String id = um.get(sym + "-" + CollectionUtil.merge(sites, "-"));

				double signedP = Double.parseDouble(t[pInd]);
				if (t[fcInd].startsWith("0")) signedP = -signedP;

				sb.append("\n").append(id).append("\t").append(sym).append("\t").append(CollectionUtil.merge(sites, "|")).append("\t").append("P").append("\t\t").append(signedP);
			}
			else cnt[mappedSites.size() > 1 ? ONE_TO_MANY : ONE_TO_NONE]++;
		});

		System.out.println("cnt[ONE_TO_ONE] = " + cnt[ONE_TO_ONE]);
		System.out.println("cnt[ONE_TO_MANY] = " + cnt[ONE_TO_MANY]);
		System.out.println("cnt[ONE_TO_NONE] = " + cnt[ONE_TO_NONE]);
		System.out.println("cnt[NO_HUM_SYM] = " + cnt[NO_HUM_SYM]);

		return sb.toString();
	}

	public static String convertTotProt(String inFile) {
		System.out.println("Newton.convertTotProt");
		StringBuilder sb = new StringBuilder();

		String[] header = FileUtil.readHeader(inFile);
		int symInd = ArrayUtil.indexOf(header, "Description");
		int[] valInds = new int[valColNames.length];
		for (int i = 0; i < valInds.length; i++)
		{
			valInds[i] = ArrayUtil.indexOf(header, valColNames[i]);
		}

		UniqueMaker um = new UniqueMaker();

		int[] cnt = new int[3];
		int ONE_TO_ONE = 0;
		int ONE_TO_NONE = 1;
		int ONE_TO_MANY = 2;

		FileUtil.linesTabbedSkip1(inFile).forEach(t ->
		{
			String sym = t[symInd];
			if (!sym.contains("OS=Mus musculus") || !sym.contains("GN=")) return;
			int beginIndex = sym.indexOf("GN=") + 3;
			int endIndex = sym.indexOf(" ", beginIndex);
			sym = endIndex < 0 ? sym.substring(beginIndex) : sym.substring(beginIndex, endIndex);

			Set<String> hSyms = MGI.get().getCorrespondingHumanSymbols(sym);

			if (hSyms.isEmpty())
			{
				cnt[ONE_TO_NONE]++;
			}
			else if (hSyms.size() > 1)
			{
				cnt[ONE_TO_MANY]++;
			}
			else
			{
				cnt[ONE_TO_ONE]++;

				sym = hSyms.iterator().next();
				String id = um.get(sym);

				sb.append("\n").append(id).append("\t").append(sym).append("\t\tG\t");
				for (int ind : valInds)
				{
					sb.append("\t").append(Math.log(Double.parseDouble(t[ind])) / Math.log(2));
				}
			}

		});
		System.out.println("cnt[ONE_TO_ONE] = " + cnt[ONE_TO_ONE]);
		System.out.println("cnt[ONE_TO_MANY] = " + cnt[ONE_TO_MANY]);
		System.out.println("cnt[ONE_TO_NONE] = " + cnt[ONE_TO_NONE]);
		return sb.toString();
	}

}
