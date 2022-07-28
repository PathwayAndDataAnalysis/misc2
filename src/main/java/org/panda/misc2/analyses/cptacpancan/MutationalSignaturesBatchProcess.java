package org.panda.misc2.analyses.cptacpancan;

import com.github.jsonldjava.utils.JsonUtils;
import org.panda.causalpath.network.GraphWriter;
import org.panda.causalpath.run.JasonizeResultGraphsRecursively;
import org.panda.misc2.causalpath.CausalPathSubnetwork;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.StreamDirection;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class MutationalSignaturesBatchProcess
{
	public static void main(String[] args) throws IOException
	{
//		processRecursive(
//			"/home/ozgunbabur/Data/CPTAC-PanCan/mutational-signatures",
//			"/home/ozgunbabur/Analyses/CPTAC-PanCan/mutational-signatures");

		addSubgraphs("/home/ozgunbabur/Analyses/CPTAC-PanCan/mutational-signatures", readGOI());
//		addSubgraphs("/home/ozgunbabur/Analyses/CPTAC-PanCan/mutational-signatures-with-kl-input", readGOI());
		jsonizeSubgraphs();

//		generateASpecificNeighborhoodGraph();
//		transposeMutexMatrix();
	}

	private static void processRecursive(String inRoot, String outRoot) throws IOException
	{
		FileUtil.processDirsRecursive(new File(inRoot), dir ->
		{
			File deFile = findTheDEFile(dir);
			if (deFile != null)
			{
				separateRegularAndResTaggedPTMData(deFile.getPath());

				CPTACPanCanSingleDEFileAnalysis.process(dir.getPath() + File.separator + "difexp.tsv", inRoot, outRoot);
				CPTACPanCanSingleDEFileAnalysis.process(dir.getPath() + File.separator + "difexp_res.tsv", inRoot, outRoot);
			}
		});
	}

	private static File findTheDEFile(File dir)
	{
		for (File file : dir.listFiles())
		{
			if (file.getName().contains("diffexp_results") && file.getName().endsWith(".tsv")) return file;
		}
		return null;
	}
	private static void separateRegularAndResTaggedPTMData(String file)
	{
		String dir = file.substring(0, file.lastIndexOf(File.separator)+1);

		String outFile1 = dir + "difexp.tsv";
		String outFile2 = dir + "difexp_res.tsv";
		BufferedWriter writer1 = FileUtil.newBufferedWriter(outFile1);
		BufferedWriter writer2 = FileUtil.newBufferedWriter(outFile2);
		String[] header = FileUtil.readHeader(file);
		int ind = ArrayUtil.indexOf(header, "feature");

		FileUtil.lines(file).forEach(l ->
		{
			if (l.startsWith("\tgene_name")) l = "index" + l;
			String[] t = l.split("\t");

			if (t[ind].equals("phosphoproteome") || t[ind].equals("acetylome"))
			{
				FileUtil.writeln(l, writer1);
			}
			else if (t[ind].equals("phosphoproteome_res") || t[ind].equals("acetylome_res"))
			{
				FileUtil.writeln(l, writer2);
			}
			else
			{
				FileUtil.writeln(l, writer1);
				FileUtil.writeln(l, writer2);
			}
		});

		FileUtil.closeWriter(writer1);
		FileUtil.closeWriter(writer2);
	}

	private static Map<String, Set<String>> readGOI() throws IOException
	{
		Map map = (Map) JsonUtils.fromInputStream(new FileInputStream(
			"/home/ozgunbabur/Data/CPTAC-PanCan/full_geneset_v3.json"));

		Map<String, Set<String>> result = new HashMap<>();

		map.forEach((o, o2) -> result.put(o.toString(), new HashSet<>((Collection) o2)));

		return result;
	}

	private static void addSubgraphs(String root, Map<String, Set<String>> subMap) throws IOException
	{
		FileUtil.processDirsRecursive(new File(root), dir ->
		{
			if (FileUtil.exists(dir.getPath() + File.separator + "causative.sif"))
			{
				subMap.forEach((name, genes) ->
				{
					try {
						CausalPathSubnetwork.writeGOINeighForCompBased(dir.getPath(), genes, StreamDirection.BOTHSTREAM, name);
					} catch (IOException e) {e.printStackTrace();}
				});
			}
		});
	}

	private static void jsonizeSubgraphs() throws IOException
	{
		Set<String> subs = new HashSet<>(readGOI().keySet());
		subs.add("causative");
		String inBase = "/home/ozgunbabur/Analyses/CPTAC-PanCan";
		JasonizeResultGraphsRecursively.generate(inBase, inBase + "/mutational-signatures", subs, inBase + "/graphs", "causative.json");
//		JasonizeResultGraphsRecursively.generate(inBase, inBase + "/mutational-signatures-with-kl-input", subs, inBase + "/graphs", "causative.json");
	}

	private static void generateASpecificNeighborhoodGraph() throws IOException
	{
		String gene = "PRKDC";
		String dir = "/home/ozgunbabur/Analyses/CPTAC-PanCan/mutational-signatures/HRD/Dendrogroup3Bvs2A/Full_GeneSpace/diff_expr_res/difexp/Group3Bvs2A_Group3B/all";
		String outSIFNoExt = gene + "-neighborhood";
		CausalPathSubnetwork.writeGOINeighForCompBased(dir, Collections.singleton(gene), StreamDirection.BOTHSTREAM, outSIFNoExt);
		String jsonDir = dir.replace("/mutational-signatures/", "/graphs/mutational-signatures/") + "/" + outSIFNoExt;
		Files.createDirectories(Paths.get(jsonDir));
		GraphWriter.convertSIFToJSON(dir + "/" + outSIFNoExt + ".sif", dir + "/" + outSIFNoExt + ".format", jsonDir + File.separator + "causative.json");
	}

	// Temporary code
	private static void transposeMutexMatrix()
	{
		String inFile = "/home/ozgunbabur/Data/CPTAC-PanCan/mutex/mmrd_ddr_matrix_binary.tsv";
		String outFile = "/home/ozgunbabur/Analyses/CPTAC-PanCan/mutex/DDR/data.txt";

		FileUtil.transpose(inFile, "\t", outFile, "\t", null, null);
	}
}
