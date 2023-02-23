package org.panda.misc2.analyses;

import org.panda.resource.HGNC;
import org.panda.resource.MGI;
import org.panda.resource.SiteMappingMouseToHuman;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.UniqueMaker;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class TanEtAL_TCR
{
	public static final String DIR = "/home/ozgunbabur/Analyses/TanEtAl_TCR/";
	public static final String[] time = new String[]{"log2(2h/0h)", "log2(8h/0h)", "log2(16h/0h)"};

	public static void main(String[] args) throws IOException
	{
		convert();
	}

	public static void convert() throws IOException
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(DIR + "CP/data.tsv");
		writer.write("ID\tSymbols\tSites\tFeature\tEffect");
		for (int i = 0; i < time.length; i++)
		{
			writer.write("\t" + time[i]);
		}

		writer.write(convertPhospho(DIR + "data/phosphoproteome.tsv"));
		writer.write(convertTotProt(DIR + "data/proteome.tsv"));
		writer.close();
	}

	public static String convertPhospho(String inFile) {
		StringBuilder sb = new StringBuilder();

		String[] header = FileUtil.readHeader(inFile, 2);
		int upInd = ArrayUtil.indexOf(header, "Protein Accession");
		int siteInd = ArrayUtil.indexOf(header, "Phosphosites");

		int[] timeInd = new int[time.length];
		for (int i = 0; i < time.length; i++)
		{
			timeInd[i] = ArrayUtil.indexOf(header, time[i]);
		}

		UniqueMaker um = new UniqueMaker();

		int[] cnt = new int[4];
		int ONE_TO_ONE = 0;
		int ONE_TO_NONE = 1;
		int ONE_TO_MANY = 2;
		int NO_HUM_SYM = 3;

		FileUtil.linesTabbedSkip(inFile, 4).forEach(t ->
		{
			String mup = t[upInd].split("\\|")[1];


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

				sb.append("\n").append(id).append("\t").append(sym).append("\t").append(CollectionUtil.merge(sites, "|")).append("\t").append("P").append("\t");
				for (int i = 0; i < time.length; i++)
				{
					sb.append("\t").append(t[timeInd[i]]);
				}
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
		StringBuilder sb = new StringBuilder();

		String[] header = FileUtil.readHeader(inFile, 2);
		int symInd = ArrayUtil.indexOf(header, "Symbol");
		int[] timeInd = new int[time.length];
		for (int i = 0; i < timeInd.length; i++)
		{
			timeInd[i] = ArrayUtil.indexOf(header, time[i]);
		}

		UniqueMaker um = new UniqueMaker();

		int[] cnt = new int[3];
		int ONE_TO_ONE = 0;
		int ONE_TO_NONE = 1;
		int ONE_TO_MANY = 2;

		FileUtil.linesTabbedSkip(inFile, 4).forEach(t ->
		{
			String sym = t[symInd];

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
				for (int ind : timeInd)
				{
					sb.append("\t").append(t[ind]);
				}
			}

		});
		System.out.println("cnt[ONE_TO_ONE] = " + cnt[ONE_TO_ONE]);
		System.out.println("cnt[ONE_TO_MANY] = " + cnt[ONE_TO_MANY]);
		System.out.println("cnt[ONE_TO_NONE] = " + cnt[ONE_TO_NONE]);
		return sb.toString();
	}

}
