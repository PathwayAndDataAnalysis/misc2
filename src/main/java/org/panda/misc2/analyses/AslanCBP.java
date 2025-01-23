package org.panda.misc2.analyses;

import org.panda.resource.HGNC;
import org.panda.resource.UniProtSequence;
import org.panda.resource.network.UniProt;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.util.HashSet;
import java.util.Set;

public class AslanCBP
{
	public static final String DATA_DIR = "/home/ozgunbabur/Data/Aslan/CBP/";
	public static final String ANALYSIS_DIR = DATA_DIR.replaceFirst("Data", "Analyses");
	public static final String ORIG_FILE = DATA_DIR + "Aslan_CBP.tsv";
	public static final String CP_DATA = ORIG_FILE.replaceFirst("Data", "Analyses");//.replaceFirst("_CBP.tsv", "_CBP_fdr.tsv");

	public static void main(String[] args)
	{
		parseData();
	}

	public static void parseData()
	{
		String[] header = FileUtil.readHeader(ORIG_FILE, 4);
		int uPInd = ArrayUtil.indexOf(header, "Protein Accessions");
		int siteInd = ArrayUtil.indexOf(header, "Site List");
		int pInd = ArrayUtil.indexOf(header, "PValue_Bypass");
//		int pInd = ArrayUtil.indexOf(header, "FDR_Bypass");
		int fcInd = ArrayUtil.indexOf(header, "FC_Bypass");

		Set<String> noMatch = new HashSet<>();

		FileUtil.createDirectories(ANALYSIS_DIR);
		BufferedWriter writer = FileUtil.newBufferedWriter(CP_DATA);
		FileUtil.write("ID\tSymbols\tSites\tEffect\tFeature\tSignedP", writer);

		FileUtil.linesTabbedSkip(ORIG_FILE, 5).forEach(t ->
		{
			Set<String> sites = new HashSet<>();
			for (String site : t[siteInd].split("; "))
			{
				sites.add(site);
			}

			Set<String> syms = new HashSet<>();
			for (String upID : t[uPInd].split("; "))
			{
				if (sitesMatch(upID, sites))
				{
					String sym = HGNC.get().getSymbol(upID);
					if (sym != null) syms.add(sym);
				}
			}
			if (syms.isEmpty())
			{
				noMatch.add(t[uPInd]);
			}
			else
			{
				double p = Double.parseDouble(t[pInd]);
				if (p == 0) p = 1E-10;
				int sign = t[fcInd].startsWith("-") ? -1 : 1;
				double signedP = sign * p;

				String symbolsStr = CollectionUtil.merge(syms, " ");
				String oneSites = CollectionUtil.merge(sites, "|");
				String siteStr = oneSites;
				for (int i = 1; i < syms.size(); i++)
				{
					siteStr += " " + oneSites;
				}
				String id = CollectionUtil.merge(syms, "_") + "_" + oneSites.replaceAll("\\|", "_");
				FileUtil.lnwrite(id + "\t" + symbolsStr + "\t" + siteStr + "\t\tP\t" + signedP, writer);
			}
		});
		FileUtil.closeWriter(writer);
	}

	static boolean sitesMatch(String upID, Set<String> sites)
	{
		for (String siteStr : sites)
		{
			int loc = Integer.parseInt(siteStr.substring(1));
			String aaa = siteStr.substring(0, 1);
			String aa = UniProtSequence.get().getAminoacidAt(upID, loc);

			if (aa == null || !aa.equals(aaa))
			{
				return false;
			}
		}
		return true;
	}
}
