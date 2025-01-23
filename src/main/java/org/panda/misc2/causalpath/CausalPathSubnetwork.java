package org.panda.misc2.causalpath;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.resource.network.PathwayCommons;
import org.panda.utility.*;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.UndirectedGraph;
import org.panda.utility.statistics.FDR;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Imagine you have a comparison-based CausalPath result with a lot of relations. It is not easily readable. So we
 * visualize a subset of those relations. Since some of those subset graphs are so small, we can fit in further PPI
 * relations to the visualizations with the same proteomic data overlay on the interacting proteins. The result is a
 * mixture of CausalPath graph and PPI networks, with all nodes showing significant proteomic data on them.
 */
public class CausalPathSubnetwork
{
	public static void main(String[] args) throws IOException
	{
//		Set<String> goi = new HashSet<>(Arrays.asList("MAPK1", "MAPK3"));
//		String dir = "/home/ozgunbabur/Analyses/Aslan-Thrombin-PAR/Regular-CP/strict-sitematch/Resting-vs-Thrombin";
//		writeGOINeighForCompBased(dir, goi, StreamDirection.BOTHSTREAM, "MAPK1-3-neigh");

//		writeChemSubset("/home/ozgunbabur/Analyses/Platelet-Blood-paper/ghatge-metabolome/cond2-relax1aa/", "causative", "causative-chem-subset");
		writeDownstreamSubgraph("/home/ozgunbabur/Analyses/Platelet-Blood-paper/ghatge-metabolome/cond2-relax1aa/", Collections.singleton("BTK"), "BTK-downstream");
	}

	public static UndirectedGraph loadPPIGraph()
	{
		UndirectedGraph graph = (UndirectedGraph) PathwayCommons.get().getGraph(SIFEnum.INTERACTS_WITH);
		graph.merge((UndirectedGraph) PathwayCommons.get().getGraph(SIFEnum.IN_COMPLEX_WITH));
		return graph;
	}

	public static void generateSubsetWithPPI(String sifFile, Collection<String> seed, UndirectedGraph ppiGraph)
		throws IOException
	{
		StringBuilder sb = new StringBuilder(sifFile.substring(0, sifFile.lastIndexOf("."))).append("-subgraph");
		seed.stream().sorted().forEach(gene -> sb.append("-").append(gene));
		sb.append(".sif");

		SIFFileUtil.generateSubsetAndAddPPI(sifFile, seed, sb.toString(), ppiGraph);
	}



	public static Set<String>[] getSignificantGenes(String sigFile, double fdrThr) throws IOException
	{
		Map<String, Double> actMap = new HashMap<>();
		Map<String, Double> inhMap = new HashMap<>();

		Files.lines(Paths.get(sigFile)).skip(2).map(l -> l.split("\t")).forEach(t ->
		{
			actMap.put(t[0], Double.valueOf(t[2]));
			inhMap.put(t[0], Double.valueOf(t[3]));
		});

		return new Set[]{new HashSet<>(FDR.select(actMap, null, fdrThr)), new HashSet<>(FDR.select(inhMap, null, fdrThr))};
	}

	public static Set<String> combine(Set<String>[] sets)
	{
		Set<String> comb = new HashSet<>(sets[0]);
		comb.addAll(sets[1]);
		return comb;
	}

	public static void generateNeighborhoodSubgraphsForSignificantsRecursively(String dir, double netSig) throws IOException
	{
		String sifFile = dir + "/causative.sif";
		String sigFile = dir + "/significance-pvals.txt";
		if (Files.exists(Paths.get(sifFile)) && Files.exists(Paths.get(sigFile)))
		{
			writeSignifNeighForCompBased(dir, StreamDirection.DOWNSTREAM, netSig);
		}
		else
		{
			for (File subdir : new File(dir).listFiles())
			{
				if (subdir.isDirectory())
				{
					generateNeighborhoodSubgraphsForSignificantsRecursively(subdir.getPath(), netSig);
				}
			}
		}
	}

	public static void printSignificantGenesRecursively(String dir, double netSig) throws IOException
	{
		printSignificantGenesRecursively(dir, dir, netSig);
	}

	public static void printSignificantGenesRecursively(String dir, String root, double netSig) throws IOException
	{
		String sigFile = dir + "/significance-pvals.txt";
		if (Files.exists(Paths.get(sigFile)))
		{
			Set<String>[] genes = getSignificantGenes(sigFile, netSig);

			Set<String> comb = combine(genes);

			List<String> sorted = comb.stream().sorted().map(g -> genes[0].contains(g) && genes[1].contains(g) ? g + "(a/i)" : genes[0].contains(g) ? g + "(a)" : g + "(i)").collect(Collectors.toList());

			if (!sorted.isEmpty())
			{
				System.out.println(dir.replace(root, "") + "\t" + CollectionUtil.merge(sorted, ", "));
			}
		}
		else
		{
			for (File subdir : new File(dir).listFiles())
			{
				if (subdir.isDirectory())
				{
					printSignificantGenesRecursively(subdir.getPath(), root, netSig);
				}
			}
		}
	}

	// Comparison based results

	public static void writeSignifNeighForCompBased(String dir, StreamDirection d, double fdrThr) throws IOException
	{
		if (!Files.exists(Paths.get(dir + "/results.txt"))) return;

		System.out.println("dir = " + dir);
		Set<String> goi = combine(getSignificantGenes(dir + "/significance-pvals.txt", fdrThr));
		System.out.println("signif = " + goi);
		System.out.println("signif size = " + goi.size());
		SIFFileUtil.writeNeighborhood(dir + "/causative.sif", goi, dir + "/causative-sig-neigh.sif", d);
		Set<String> ids = getIDsAtTheNeighborhood(dir + "/results.txt", goi, d);
		System.out.println("ids.size() = " + ids.size());
		writeSubsetFormat(dir + "/causative.format", dir + "/causative-sig-neigh.format", null, ids);
	}

	public static void writeGOINeighForCompBasedRecursively(String dir, Set<String> goi, StreamDirection d) throws IOException
	{
		writeGOINeighForCompBasedRecursively(dir, goi, d, "causative-goi-neigh");
	}

	public static void writeGOINeighForCompBasedRecursively(String dir, Set<String> goi, StreamDirection d, String outSIFNoExt) throws IOException
	{
		String sifFile = dir + "/causative.sif";
		if (Files.exists(Paths.get(sifFile)))
		{
			writeGOINeighForCompBased(dir, goi, d, outSIFNoExt);
		}
		else
		{
			for (File subdir : new File(dir).listFiles())
			{
				if (subdir.isDirectory())
				{
					writeGOINeighForCompBasedRecursively(subdir.getPath(), goi, d, outSIFNoExt);
				}
			}
		}
	}

	public static void writeGOINeighForCompBasedRecursively(String dir, Map<String, Set<String>> goiMap, StreamDirection d) throws IOException
	{
		String sifFile = dir + "/causative.sif";
		if (Files.exists(Paths.get(sifFile)))
		{
			writeGOINeighForCompBased(dir, goiMap, d);
		}
		else
		{
			for (File subdir : new File(dir).listFiles())
			{
				if (subdir.isDirectory())
				{
					writeGOINeighForCompBasedRecursively(subdir.getPath(), goiMap, d);
				}
			}
		}
	}

	public static void writeGOINeighForCompBased(String dir, Set<String> goi, StreamDirection d) throws IOException
	{
		writeGOINeighForCompBased(dir, goi, d, "causative-goi-neigh");
	}

	public static void writeGOINeighForCompBased(String dir, Map<String, Set<String>> goiMap, StreamDirection d) throws IOException
	{
		for (String name : goiMap.keySet())
		{
			writeGOINeighForCompBased(dir, goiMap.get(name), d, name);
		}
	}

	public static void writeGOINeighForCompBased(String dir, Set<String> goi, StreamDirection d, String outSIFNoExt) throws IOException
	{
		if (!Files.exists(Paths.get(dir + "/results.txt"))) return;

		System.out.println("dir = " + dir);
		if (SIFFileUtil.writeNeighborhood(dir + "/causative.sif", goi, dir + "/" + outSIFNoExt + ".sif", d))
		{
			Set<String> ids = getIDsAtTheNeighborhood(dir + "/results.txt", goi, d);
			System.out.println("ids.size() = " + ids.size());
//			writeSubsetFormat(dir + "/causative.format", dir + "/" + outSIFNoExt + ".format", goi, null);
//			writeSubsetFormat(dir + "/causative.format", dir + "/" + outSIFNoExt + ".format", null, ids);
			writeSubsetFormat(dir + "/causative.format", dir + "/" + outSIFNoExt + ".format", null, null);
		}
	}

	public static void writeChemSubset(String dir, String inSIFNoExt, String outSIFNoExt) throws IOException
	{
		BufferedWriter writer = FileUtil.newBufferedWriter(dir + "/" + outSIFNoExt + ".sif");
		FileUtil.lines(dir + "/" + inSIFNoExt + ".sif").filter(l -> l.contains("CHEBI:")).forEach(l -> FileUtil.writeln(l, writer));
		writer.close();
		writeSubsetFormat(dir + "/" + inSIFNoExt + ".format", dir + "/" + outSIFNoExt + ".format", null, null);
	}

	public static void writeDownstreamSubgraph(String dir, Set<String> seeds, String outSIFNoExt) throws IOException
	{
		if (!Files.exists(Paths.get(dir + "/results.txt"))) return;

		System.out.println("dir = " + dir);

		DirectedGraph graph = SIFFileUtil.loadDirectedSingleGraph(dir + "/causative.sif");
		DirectedGraph subgraph = graph.getDownstreamSubgraph(seeds);

		if (subgraph.getSymbols().isEmpty()) return;

		SIFFileUtil.writeSubgraph(dir + "/causative.sif", subgraph, dir + "/" + outSIFNoExt + ".sif");
		Set<String> ids = getIDsInTheSubgraph(dir + "/results.txt", subgraph);
		System.out.println("ids.size() = " + ids.size());
		writeSubsetFormat(dir + "/causative.format", dir + "/" + outSIFNoExt + ".format", null, ids);
	}

	public static void writeIntersectionSubgraph(String primaryDir, String primarySIFName, String secondaryFile, String outSIFNoExt) throws IOException
	{
		if (!Files.exists(Paths.get(primaryDir + "/results.txt"))) return;

		System.out.println("dir = " + primaryDir);

		Set<String> keep = FileUtil.linesTabbed(secondaryFile).filter(t -> t.length > 2).map(t -> t[0] + "\t" + t[1] + "\t" + t[2]).collect(Collectors.toSet());


		SIFFileUtil.writeSubgraph(primaryDir + "/" + primarySIFName + ".sif", keep, primaryDir + "/" + outSIFNoExt + ".sif");
		Set<String> ids = getIDsInTheSubgraph(primaryDir + "/results.txt", keep);
		System.out.println("ids.size() = " + ids.size());
		writeSubsetFormat(primaryDir + "/" + primarySIFName + ".format", primaryDir + "/" + outSIFNoExt + ".format", null, ids);
	}

	// Correlation-based results

	/**
	 * @deprecated Provide the FDR threshold and use the next method.
	 */
	public static void writeSignifNeighForCorrBased(String dir, StreamDirection d) throws IOException
	{
		writeSignifNeighForCorrBased(dir, d, 0.1);
	}

	public static void writeSignifNeighForCorrBased(String dir, StreamDirection d, double fdrThr) throws IOException
	{
		System.out.println("dir = " + dir);
		Set<String> goi = getDownstreamEnrichedForCorrelation(dir + "/significance-pvals.txt", fdrThr);
		System.out.println("signif size = " + goi.size());
		SIFFileUtil.writeNeighborhood(dir + "/causative.sif", goi, dir + "/causative-sig-neigh.sif", d);
		Set<String> ids = getIDsAtTheNeighborhood(dir + "/results.txt", goi, d);
		System.out.println("ids.size() = " + ids.size());
		writeSubsetFormat(dir + "/causative.format", dir + "/causative-sig-neigh.format", ids);
	}

	public static void writeSignifNeighForCorrBasedRecursive(String base, StreamDirection d, double fdrThr) throws IOException
	{
		FileUtil.processDirsRecursive(new File(base), dir ->
		{
			if (Files.exists(Paths.get(dir.getPath() + "/significance-pvals.txt")))
			{
				writeSignifNeighForCorrBased(dir.getPath(), d, fdrThr);
			}
		});
	}

	public static void writeGOINeighForCorrBased(String dir, Set<String> goi, StreamDirection d) throws IOException
	{
		writeGOINeighForCorrBased(dir, goi, d, "goi-" + d.toString().toLowerCase());
	}

	public static void writeGOINeighForCorrBased(String dir, Set<String> goi, StreamDirection d, String resultSIFNoExt) throws IOException
	{
		writeGOINeighForCorrBased(dir, goi, d, resultSIFNoExt, true);
	}

	public static void writeGOINeighForCorrBased(String dir, Set<String> goi, StreamDirection d, String resultSIFNoExt, boolean markGOI) throws IOException
	{
		System.out.println("dir = " + dir);
		SIFFileUtil.writeNeighborhood(dir + "/causative.sif", goi, dir + "/" +resultSIFNoExt + ".sif", d);
		Set<String> ids = getIDsAtTheNeighborhood(dir + "/results.txt", goi, d);
		System.out.println("ids.size() = " + ids.size());
		writeSubsetFormat(dir + "/causative.format", dir + "/" + resultSIFNoExt + ".format", markGOI ? goi : null, ids);
	}

	public static void writeSubsetFormat(String inFile, String outFile, Set<String> markSubset, Set<String> ids) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		Files.lines(Paths.get(inFile)).forEach(l ->
		{
			String[] t = l.split("\t");

			if (t[2].equals("rppasite"))
			{
				if (ids == null || ids.contains(t[3].split("\\|")[0]))
				{
					FileUtil.writeln(l, writer);
				}
			}
			else FileUtil.writeln(l, writer);
		});

		if (markSubset != null) markSubset.forEach(g -> FileUtil.writeln("node\t" + g + "\tbordercolor\t100 0 255", writer));

		writer.close();
	}

	private static void writeSubsetFormat(String inFile, String outFile, Set<String> ids) throws IOException
	{
		writeSubsetFormat(inFile, outFile, null, ids);
	}

	private static Set<String> getDownstreamEnrichedForCorrelation(String sigFile, double fdrThr) throws IOException
	{
		Map<String, Double> pMap = Files.lines(Paths.get(sigFile)).skip(2).map(l -> l.split("\t")).collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));
		return new HashSet<>(FDR.select(pMap, null, fdrThr));
	}

	public static Set<String> getIDsAtTheNeighborhood(String resultFile, Set<String> goi, StreamDirection d) throws IOException
	{
		String[] header = Files.lines(Paths.get(resultFile)).findFirst().get().split("\t");
		int sInd = ArrayUtil.indexOf(header, "Source");
		int tInd = ArrayUtil.indexOf(header, "Target");
		int sDI = ArrayUtil.indexOf(header, "Source data ID");
		int tDI = ArrayUtil.indexOf(header, "Target data ID");

		Set<String> ids = new HashSet<>();

		Files.lines(Paths.get(resultFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if ((d == StreamDirection.DOWNSTREAM && goi.contains(t[sInd])) ||
				(d == StreamDirection.UPSTREAM && goi.contains(t[tInd])) ||
				(d == StreamDirection.BOTHSTREAM && (goi.contains(t[tInd]) || goi.contains(t[sInd]))))
			{
				ids.add(t[sDI]);
				ids.add(t[tDI]);
			}
		});
		return ids;
	}

	public static Set<String> getIDsInTheSubgraph(String resultFile, DirectedGraph subgraph) throws IOException
	{
		String[] header = Files.lines(Paths.get(resultFile)).findFirst().get().split("\t");
		int sInd = ArrayUtil.indexOf(header, "Source");
		int tInd = ArrayUtil.indexOf(header, "Target");
		int sDI = ArrayUtil.indexOf(header, "Source data ID");
		int tDI = ArrayUtil.indexOf(header, "Target data ID");

		Set<String> ids = new HashSet<>();

		Files.lines(Paths.get(resultFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if (subgraph.hasRelation(t[sInd], t[tInd]))
			{
				ids.add(t[sDI]);
				ids.add(t[tDI]);
			}
		});
		return ids;
	}

	/**
	 * the select set has to contain strings that have source, relation and target as tab separated.
	 */
	public static Set<String> getIDsInTheSubgraph(String resultFile, Set<String> select) throws IOException
	{
		String[] header = Files.lines(Paths.get(resultFile)).findFirst().get().split("\t");
		int sInd = ArrayUtil.indexOf(header, "Source");
		int tInd = ArrayUtil.indexOf(header, "Target");
		int rInd = ArrayUtil.indexOf(header, "Relation");
		int sDI = ArrayUtil.indexOf(header, "Source data ID");
		int tDI = ArrayUtil.indexOf(header, "Target data ID");

		Set<String> ids = new HashSet<>();

		Files.lines(Paths.get(resultFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			if (select.contains(t[sInd] + "\t" + t[rInd] + "\t" + t[tInd]))
			{
				ids.add(t[sDI]);
				ids.add(t[tDI]);
			}
		});
		return ids;
	}

	/**
	 * It is assumed that the goiDir directory contains text files that have gene symbols in them (one per line).
	 * @param parentDir
	 * @param goiDir
	 */
	public static void writeGOIRecursiveCompBased(String parentDir, String goiDir) throws IOException
	{
		Map<String, Set<String>> goiMap = new HashMap<>();
		for (File file : new File(goiDir).listFiles())
		{
			String name = file.getName();
			if (name.endsWith(".txt"))
			{
				Set<String> goi = Files.lines(Paths.get(file.getPath())).collect(Collectors.toSet());
				goiMap.put(name.substring(0, name.length() - 4), goi);
			}
		}

		writeGOIRecursiveCompBased(parentDir, goiMap);
	}

	public static void writeGOIRecursiveCompBased(String dir, Map<String, Set<String>> goiMap) throws IOException
	{
		if (Files.exists(Paths.get(dir + "/results.txt")))
		{
			System.out.println("dir = " + dir);
			for (String goiName : goiMap.keySet())
			{
				Set<String> goi = goiMap.get(goiName);
				writeGOINeighForCompBased(dir, goi, StreamDirection.BOTHSTREAM, "causative-" + goiName);
			}
		}

		for (File child : new File(dir).listFiles())
		{
			if (child.isDirectory()) writeGOIRecursiveCompBased(child.getPath(), goiMap);
		}
	}

	public static void writeSiteSpecAndStrongExpCompBased(String dir, String outWoExt) throws IOException
	{
		Map<String, Integer> rnaDir = FileUtil.linesTabbed(dir + "/results.txt").skip(1)
			.filter(t -> t[7].endsWith("-rna")).collect(Collectors.toMap(
				t -> t[2], t -> t[8].startsWith("-") ? -1 : 1, (i1, i2) -> i1));

		Map<String, Integer> prtDir = FileUtil.linesTabbed(dir + "/results.txt").skip(1)
			.filter(t -> t[7].equals(t[2])).collect(Collectors.toMap(
				t -> t[2], t -> t[8].startsWith("-") ? -1 : 1, (i1, i2) -> i1));

		Set<String> consistent = rnaDir.keySet().stream().filter(prtDir.keySet()::contains)
			.filter(g -> rnaDir.get(g).equals(prtDir.get(g))).collect(Collectors.toSet());

		Set<String> sigGenes = FileUtil.linesTabbed(dir + "/results.txt").skip(1)
			.filter(t -> t[4].endsWith("active-by-network-sig")).map(t -> t[0]).collect(Collectors.toSet());

		BufferedWriter writer1 = FileUtil.newBufferedWriter(dir + "/" + outWoExt + ".sif");
		FileUtil.linesTabbed(dir + "/causative.sif").filter(t -> t.length > 2 && t[1].contains("expression")
			&& (sigGenes.contains(t[0]) && consistent.contains(t[2])))
			.forEach(t -> FileUtil.writeln(ArrayUtil.merge("\t", t), writer1));
		writer1.close();

		FileUtil.copyFile(dir + "/causative.format", dir + "/" + outWoExt + ".format");
	}
}
