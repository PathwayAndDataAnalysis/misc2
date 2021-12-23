package org.panda.misc2.analyses.cptacpancan;

import org.panda.causalpath.run.JasonizeResultGraphsRecursively;
import org.panda.misc2.causalpath.CausalPathSubnetwork;
import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class KinaseLibraryIntegrator
{
	public static void main(String[] args) throws IOException
	{
//		prepareCPDirs();
		compareCPWithKL();
//		jsonizeSubgraphs();
	}

	private static void prepareCPDirs() throws IOException
	{
		String inBase = "/home/ozgunbabur/Analyses/CPTAC-PanCan/mutational-signatures/";
		String outBase = "/home/ozgunbabur/Analyses/CPTAC-PanCan/mutational-signatures-with-kl-input/";
		String klBase = "/home/ozgunbabur/Data/CPTAC-PanCan/KinaseLibrary/";

		// HRD vs HRP

		prepareCPDir(inBase + "HRD/HRDvsHRP/DDR_subset/diff_expr_res/difexp/HRD/all",
			outBase + "HRD/HRDvsHRP/DDR_subset/difexp/HRD",
			klBase + "HRD_MMRD/phosphoproteome/ser_thr/tables/HRD_vs_HRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome/tyrosine/tables/HRD_vs_HRP_FC05.tsv");

		prepareCPDir(inBase + "HRD/HRDvsHRP/DDR_subset/diff_expr_res/difexp_res/HRD/all",
			outBase + "HRD/HRDvsHRP/DDR_subset/difexp_res/HRD",
			klBase + "HRD_MMRD/phosphoproteome_res/ser_thr/tables/HRD_vs_HRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome_res/tyrosine/tables/HRD_vs_HRP_FC05.tsv");

		prepareCPDir(inBase + "HRD/HRDvsHRP/Full_GeneSpace/diff_expr_res/difexp/HRD/all",
			outBase + "HRD/HRDvsHRP/Full_GeneSpace/difexp/HRD",
			klBase + "HRD_MMRD/phosphoproteome/ser_thr/tables/HRD_vs_HRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome/tyrosine/tables/HRD_vs_HRP_FC05.tsv");

		prepareCPDir(inBase + "HRD/HRDvsHRP/Full_GeneSpace/diff_expr_res/difexp_res/HRD/all",
			outBase + "HRD/HRDvsHRP/Full_GeneSpace/difexp_res/HRD",
			klBase + "HRD_MMRD/phosphoproteome_res/ser_thr/tables/HRD_vs_HRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome_res/tyrosine/tables/HRD_vs_HRP_FC05.tsv");

		// MMRD vs MMRP

		prepareCPDir(inBase + "MMRD/MMRDvsMMRP/DDR_Subset/diff_expr_res/difexp/MMRD/all",
			outBase + "MMRD/MMRDvsMMRP/DDR_subset/difexp/MMRD",
			klBase + "HRD_MMRD/phosphoproteome/ser_thr/tables/MMRD_vs_MMRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome/tyrosine/tables/MMRD_vs_MMRP_FC05.tsv");

		prepareCPDir(inBase + "MMRD/MMRDvsMMRP/DDR_Subset/diff_expr_res/difexp_res/MMRD/all",
			outBase + "MMRD/MMRDvsMMRP/DDR_subset/difexp_res/MMRD",
			klBase + "HRD_MMRD/phosphoproteome_res/ser_thr/tables/MMRD_vs_MMRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome_res/tyrosine/tables/MMRD_vs_MMRP_FC05.tsv");

		prepareCPDir(inBase + "MMRD/MMRDvsMMRP/Full_GeneSet/diff_expr_res/difexp/MMRD/all",
			outBase + "MMRD/MMRDvsMMRP/Full_GeneSet/difexp/MMRD",
			klBase + "HRD_MMRD/phosphoproteome/ser_thr/tables/MMRD_vs_MMRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome/tyrosine/tables/MMRD_vs_MMRP_FC05.tsv");

		prepareCPDir(inBase + "MMRD/MMRDvsMMRP/Full_GeneSet/diff_expr_res/difexp_res/MMRD/all",
			outBase + "MMRD/MMRDvsMMRP/Full_GeneSet/difexp_res/MMRD",
			klBase + "HRD_MMRD/phosphoproteome_res/ser_thr/tables/MMRD_vs_MMRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome_res/tyrosine/tables/MMRD_vs_MMRP_FC05.tsv");

		// 3b vs 2a

		prepareCPDir(inBase + "HRD/Dendrogroup3Bvs2A/DDR_Subset/diff_expr_res/difexp/Group3Bvs2A_Group3B/all",
			outBase + "HRD/Dendrogroup3Bvs2A/DDR_Subset/difexp/Group3Bvs2A_Group3B",
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_3b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_3b_vs_2a_FC05.tsv");

		prepareCPDir(inBase + "HRD/Dendrogroup3Bvs2A/DDR_Subset/diff_expr_res/difexp_res/Group3Bvs2A_Group3B/all",
			outBase + "HRD/Dendrogroup3Bvs2A/DDR_Subset/difexp_res/Group3Bvs2A_Group3B",
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_3b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_3b_vs_2a_FC05.tsv");

		prepareCPDir(inBase + "HRD/Dendrogroup3Bvs2A/Full_GeneSpace/diff_expr_res/difexp/Group3Bvs2A_Group3B/all",
			outBase + "HRD/Dendrogroup3Bvs2A/Full_GeneSpace/difexp/Group3Bvs2A_Group3B",
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_3b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_3b_vs_2a_FC05.tsv");

		prepareCPDir(inBase + "HRD/Dendrogroup3Bvs2A/Full_GeneSpace/diff_expr_res/difexp_res/Group3Bvs2A_Group3B/all",
			outBase + "HRD/Dendrogroup3Bvs2A/Full_GeneSpace/difexp_res/Group3Bvs2A_Group3B",
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_3b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_3b_vs_2a_FC05.tsv");

		// 5b vs 2a

		prepareCPDir(inBase + "HRD/Dendrogroup5Bvs2A/DDR_Subset/difexp/Group5Bvs2A_Group5B/all",
			outBase + "HRD/Dendrogroup5Bvs2A/DDR_Subset/difexp/Group5Bvs2A_Group5B",
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_5b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_5b_vs_2a_FC05.tsv");

		prepareCPDir(inBase + "HRD/Dendrogroup5Bvs2A/DDR_Subset/difexp_res/Group5Bvs2A_Group5B/all",
			outBase + "HRD/Dendrogroup5Bvs2A/DDR_Subset/difexp_res/Group5Bvs2A_Group5B",
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_5b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_5b_vs_2a_FC05.tsv");

		prepareCPDir(inBase + "HRD/Dendrogroup5Bvs2A/Full_GeneSpace/diff_expr_res/difexp/Group5Bvs2A_Group5B/all",
			outBase + "HRD/Dendrogroup5Bvs2A/Full_GeneSpace/difexp/Group5Bvs2A_Group5B",
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_5b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_5b_vs_2a_FC05.tsv");

		prepareCPDir(inBase + "HRD/Dendrogroup5Bvs2A/Full_GeneSpace/diff_expr_res/difexp_res/Group5Bvs2A_Group5B/all",
			outBase + "HRD/Dendrogroup5Bvs2A/Full_GeneSpace/difexp_res/Group5Bvs2A_Group5B",
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_5b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_5b_vs_2a_FC05.tsv");

		// 5b vs 10b

		prepareCPDir(inBase + "HRD/Dendrogroup5Bvs10B/DDR_Subset/diff_expr_res/difexp/Group5Bvs10B_Group5B/all",
			outBase + "HRD/Dendrogroup5Bvs10B/DDR_Subset/difexp/Group5Bvs10B_Group5B",
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_5b_vs_10b_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_5b_vs_10b_FC05.tsv");

		prepareCPDir(inBase + "HRD/Dendrogroup5Bvs10B/DDR_Subset/diff_expr_res/difexp_res/Group5Bvs10B_Group5B/all",
			outBase + "HRD/Dendrogroup5Bvs10B/DDR_Subset/difexp_res/Group5Bvs10B_Group5B",
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_5b_vs_10b_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_5b_vs_10b_FC05.tsv");

		prepareCPDir(inBase + "HRD/Dendrogroup5Bvs10B/Full_GeneSpace/diff_expr_res/difexp/Group5Bvs10B_Group5B/all",
			outBase + "HRD/Dendrogroup5Bvs10B/Full_GeneSpace/difexp/Group5Bvs10B_Group5B",
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_5b_vs_10b_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_5b_vs_10b_FC05.tsv");

		prepareCPDir(inBase + "HRD/Dendrogroup5Bvs10B/Full_GeneSpace/diff_expr_res/difexp_res/Group5Bvs10B_Group5B/all",
			outBase + "HRD/Dendrogroup5Bvs10B/Full_GeneSpace/difexp_res/Group5Bvs10B_Group5B",
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_5b_vs_10b_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_5b_vs_10b_FC05.tsv");
	}

	private static void prepareCPDir(String existingCPDir, String newCPDir, String... klFiles) throws IOException
	{
		Set<String>[] sig = getSignificant(0.1, klFiles);

		if (sig[0].isEmpty() && sig[1].isEmpty()) return;

		FileUtil.mkdirs(newCPDir);
		BufferedWriter writer = FileUtil.newBufferedWriter(newCPDir + "/parameters.txt");

		String oneLess = existingCPDir.substring(0, existingCPDir.lastIndexOf(File.separator));
		String dataAbsPath = oneLess.substring(0, oneLess.lastIndexOf(File.separator)) + "/data.txt";

		FileUtil.lines(existingCPDir + "/parameters.txt")
			.map(l -> l.replaceAll("../../data.txt", dataAbsPath))
			.forEach(l -> FileUtil.writeln(l, writer));

		sig[0].forEach(g -> FileUtil.lnwrite("gene-activity = " + g + " a", writer));
		sig[1].forEach(g -> FileUtil.lnwrite("gene-activity = " + g + " i", writer));

		Set<String> sigs = new HashSet<>(sig[0]);
		sigs.addAll(sig[1]);
		FileUtil.lnwrite("gene-focus = " + CollectionUtil.merge(sigs, ";"), writer);

		writer.close();
	}


	private static Set<String>[] getSignificant(double thr, String... file)
	{
		Set<String> ups = new HashSet<>();
		Set<String> dws = new HashSet<>();

		for (String filename : file)
		{
			String[] header = FileUtil.readHeader(filename);
			int ffuIndex = ArrayUtil.indexOf(header, "freq_factor_upreg");
			int ffdIndex = ArrayUtil.indexOf(header, "freq_factor_downreg");
			int puIndex = ArrayUtil.indexOf(header, "adj.pval_upreg");
			int pdIndex = ArrayUtil.indexOf(header, "adj.pval_downreg");

			FileUtil.linesTabbedSkip1(filename).forEach(t ->
			{
				if (!t[ffuIndex].isEmpty() && Double.parseDouble(t[ffuIndex]) > 1)
				{
					double p = Double.parseDouble(t[puIndex]);
					if (p <= thr) ups.add(t[0]);
				}
				if (!t[ffdIndex].isEmpty() && (t[ffdIndex].equals("inf") || Double.parseDouble(t[ffdIndex]) > 1))
				{
					double p = Double.parseDouble(t[pdIndex]);
					if (p <= thr) dws.add(t[0]);
				}
			});
		}

		return new Set[]{ups, dws};
	}

	private static void compareCPWithKL() throws IOException
	{
		BufferedWriter writer = FileUtil.newBufferedWriter("/home/ozgunbabur/Analyses/CPTAC-PanCan/cp-kl-comparison.tsv");

		String inBase = "/home/ozgunbabur/Analyses/CPTAC-PanCan/mutational-signatures/";
		String outBase = "/home/ozgunbabur/Analyses/CPTAC-PanCan/mutational-signatures-with-kl-input/";
		String klBase = "/home/ozgunbabur/Data/CPTAC-PanCan/KinaseLibrary/";

		// HRD vs HRP

		compareCPWithKL(inBase + "HRD/HRDvsHRP/DDR_subset/diff_expr_res/difexp/HRD/all", writer,
			klBase + "HRD_MMRD/phosphoproteome/ser_thr/tables/HRD_vs_HRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome/tyrosine/tables/HRD_vs_HRP_FC05.tsv");

		compareCPWithKL(inBase + "HRD/HRDvsHRP/DDR_subset/diff_expr_res/difexp_res/HRD/all", writer,
			klBase + "HRD_MMRD/phosphoproteome_res/ser_thr/tables/HRD_vs_HRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome_res/tyrosine/tables/HRD_vs_HRP_FC05.tsv");

		compareCPWithKL(inBase + "HRD/HRDvsHRP/Full_GeneSpace/diff_expr_res/difexp/HRD/all", writer,
			klBase + "HRD_MMRD/phosphoproteome/ser_thr/tables/HRD_vs_HRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome/tyrosine/tables/HRD_vs_HRP_FC05.tsv");

		compareCPWithKL(inBase + "HRD/HRDvsHRP/Full_GeneSpace/diff_expr_res/difexp_res/HRD/all", writer,
			klBase + "HRD_MMRD/phosphoproteome_res/ser_thr/tables/HRD_vs_HRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome_res/tyrosine/tables/HRD_vs_HRP_FC05.tsv");

		// MMRD vs MMRP

		compareCPWithKL(inBase + "MMRD/MMRDvsMMRP/DDR_Subset/diff_expr_res/difexp/MMRD/all", writer,
			klBase + "HRD_MMRD/phosphoproteome/ser_thr/tables/MMRD_vs_MMRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome/tyrosine/tables/MMRD_vs_MMRP_FC05.tsv");

		compareCPWithKL(inBase + "MMRD/MMRDvsMMRP/DDR_Subset/diff_expr_res/difexp_res/MMRD/all", writer,
			klBase + "HRD_MMRD/phosphoproteome_res/ser_thr/tables/MMRD_vs_MMRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome_res/tyrosine/tables/MMRD_vs_MMRP_FC05.tsv");

		compareCPWithKL(inBase + "MMRD/MMRDvsMMRP/Full_GeneSet/diff_expr_res/difexp/MMRD/all", writer,
			klBase + "HRD_MMRD/phosphoproteome/ser_thr/tables/MMRD_vs_MMRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome/tyrosine/tables/MMRD_vs_MMRP_FC05.tsv");

		compareCPWithKL(inBase + "MMRD/MMRDvsMMRP/Full_GeneSet/diff_expr_res/difexp_res/MMRD/all", writer,
			klBase + "HRD_MMRD/phosphoproteome_res/ser_thr/tables/MMRD_vs_MMRP_FC05.tsv",
			klBase + "HRD_MMRD/phosphoproteome_res/tyrosine/tables/MMRD_vs_MMRP_FC05.tsv");

		// 3b vs 2a

		compareCPWithKL(inBase + "HRD/Dendrogroup3Bvs2A/DDR_Subset/diff_expr_res/difexp/Group3Bvs2A_Group3B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_3b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_3b_vs_2a_FC05.tsv");

		compareCPWithKL(inBase + "HRD/Dendrogroup3Bvs2A/DDR_Subset/diff_expr_res/difexp_res/Group3Bvs2A_Group3B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_3b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_3b_vs_2a_FC05.tsv");

		compareCPWithKL(inBase + "HRD/Dendrogroup3Bvs2A/Full_GeneSpace/diff_expr_res/difexp/Group3Bvs2A_Group3B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_3b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_3b_vs_2a_FC05.tsv");

		compareCPWithKL(inBase + "HRD/Dendrogroup3Bvs2A/Full_GeneSpace/diff_expr_res/difexp_res/Group3Bvs2A_Group3B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_3b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_3b_vs_2a_FC05.tsv");

		// 5b vs 2a

		compareCPWithKL(inBase + "HRD/Dendrogroup5Bvs2A/DDR_Subset/difexp/Group5Bvs2A_Group5B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_5b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_5b_vs_2a_FC05.tsv");

		compareCPWithKL(inBase + "HRD/Dendrogroup5Bvs2A/DDR_Subset/difexp_res/Group5Bvs2A_Group5B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_5b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_5b_vs_2a_FC05.tsv");

		compareCPWithKL(inBase + "HRD/Dendrogroup5Bvs2A/Full_GeneSpace/diff_expr_res/difexp/Group5Bvs2A_Group5B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_5b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_5b_vs_2a_FC05.tsv");

		compareCPWithKL(inBase + "HRD/Dendrogroup5Bvs2A/Full_GeneSpace/diff_expr_res/difexp_res/Group5Bvs2A_Group5B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_5b_vs_2a_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_5b_vs_2a_FC05.tsv");

		// 5b vs 10b

		compareCPWithKL(inBase + "HRD/Dendrogroup5Bvs10B/DDR_Subset/diff_expr_res/difexp/Group5Bvs10B_Group5B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_5b_vs_10b_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_5b_vs_10b_FC05.tsv");

		compareCPWithKL(inBase + "HRD/Dendrogroup5Bvs10B/DDR_Subset/diff_expr_res/difexp_res/Group5Bvs10B_Group5B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_5b_vs_10b_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_5b_vs_10b_FC05.tsv");

		compareCPWithKL(inBase + "HRD/Dendrogroup5Bvs10B/Full_GeneSpace/diff_expr_res/difexp/Group5Bvs10B_Group5B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome/ser_thr/tables/dendrogroup_5b_vs_10b_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome/tyrosine/tables/dendrogroup_5b_vs_10b_FC05.tsv");

		compareCPWithKL(inBase + "HRD/Dendrogroup5Bvs10B/Full_GeneSpace/diff_expr_res/difexp_res/Group5Bvs10B_Group5B/all", writer,
			klBase + "DDR-Hypoxia/phosphoproteome_res/ser_thr/tables/dendrogroup_5b_vs_10b_FC05.tsv",
			klBase + "DDR-Hypoxia/phosphoproteome_res/tyrosine/tables/dendrogroup_5b_vs_10b_FC05.tsv");

		writer.close();
	}

	private static void compareCPWithKL(String cpDir, BufferedWriter writer, String... klFiles) throws IOException
	{
		Set<String> cpGenes = FileUtil.linesTabbedSkip1(cpDir + "/causative.sif").filter(t -> t.length > 2)
			.map(t -> new HashSet<>(Arrays.asList(t[0], t[2]))).flatMap(Collection::stream).collect(Collectors.toSet());

		Set<String>[] cpSigs = CausalPathSubnetwork.getSignificantGenes(cpDir + "/significance-pvals.txt", 0.1);

		Set<String>[] klSigs = getSignificant(0.1, klFiles);

		Set<String> cpUpDif = new HashSet<>(cpSigs[0]);
		Set<String> cpDwDif = new HashSet<>(cpSigs[1]);
		Set<String> klUpDif = new HashSet<>(klSigs[0]);
		Set<String> klDwDif = new HashSet<>(klSigs[1]);
		Set<String> upCommon = new HashSet<>(cpSigs[0]);
		Set<String> dwCommon = new HashSet<>(cpSigs[1]);
		upCommon.retainAll(klSigs[0]);
		dwCommon.retainAll(klSigs[1]);
		cpUpDif.removeAll(upCommon);
		cpDwDif.removeAll(dwCommon);
		klUpDif.removeAll(upCommon);
		klDwDif.removeAll(dwCommon);

		System.out.println("upCommon = " + upCommon);
		System.out.println("dwCommon = " + dwCommon);

		FileUtil.lnwrite(cpDir + "\t" + cpUpDif.size() + "\t" + upCommon.size() + "\t" + klUpDif.size() + "\t" +
			cpDwDif.size() + "\t" + dwCommon.size() + "\t" + klDwDif.size(), writer);

		Set<String> cpAllSigs = new HashSet<>(cpSigs[0]);
		cpAllSigs.addAll(cpSigs[1]);

		Set<String> klAll = new HashSet<>(klSigs[0]);
		klAll.addAll(klSigs[1]);
		System.out.println("cpDir = " + cpDir);
		CollectionUtil.printVennSets(10, cpGenes, cpAllSigs, klAll);
	}

	private static void jsonizeSubgraphs() throws IOException
	{
		String inBase = "/home/ozgunbabur/Analyses/CPTAC-PanCan";
		JasonizeResultGraphsRecursively.generate(inBase, inBase + "/mutational-signatures-with-kl-input", Collections.singleton("causative"), inBase + "/graphs", "causative.json");
	}

}