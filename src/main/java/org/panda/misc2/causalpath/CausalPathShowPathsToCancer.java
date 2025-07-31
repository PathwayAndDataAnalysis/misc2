package org.panda.misc2.causalpath;

import org.panda.misc2.omics.IdentifyOncogenicProteomicChanges;
import org.panda.resource.OncoKB;
import org.panda.resource.signednetwork.SignedType;
import org.panda.resource.siteeffect.SiteEffectCollective;
import org.panda.resource.siteeffect.Feature;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Assumes each proteomic row has one gene.
 */
public class CausalPathShowPathsToCancer
{

	public static void main(String[] args) throws IOException
	{
//		generateRecursivelyComparison("/Users/ozgun/Documents/Analyses/CPTAC-LSCC-3.2/clusters/against-normals", Boolean.TRUE);
//		generateRecursivelyComparison("/Users/ozgun/Documents/Analyses/CPTAC-LSCC-3.2/temp", null);
		generateRecursivelyComparison("/home/ozgunbabur/Analyses/CPTAC-LSCC/v3/NMF-subtype-vs-others", false);
//		generateRecursivelyComparison("/home/ozgunbabur/Analyses/CPTAC-LSCC/v3/tumors-vs-normals-per-NMF-subtype", null);
	}

	public static void generateRecursivelyComparison(String dir, Boolean geneexp) throws IOException
	{
		if (Files.exists(Paths.get(dir + "/results.txt")))
			writeOncogenicRelations(dir, geneexp);

		for (File d : new File(dir).listFiles())
		{
			if (d.isDirectory()) generateRecursivelyComparison(d.getPath(), geneexp);
		}
	}

	public static void generateRecursivelyCorrelation(String dir) throws IOException
	{
		if (Files.exists(Paths.get(dir + "/results.txt")))
			writeOncogenicRelationsForCorrelation(dir);

		for (File d : new File(dir).listFiles())
		{
			if (d.isDirectory()) generateRecursivelyCorrelation(d.getPath());
		}
	}

	private static void writeOncogenicRelations(String dir, Boolean geneexp) throws IOException
	{
		String outName = "/oncogenic-changes" + (geneexp == null ? "" : !geneexp ? "-phospho" : "-express");
//		String outName = "/oncogenic-changes-lung" + (geneexp == null ? "" : !geneexp ? "-phospho" : "-express");

		String resultFile = dir + "/results.txt";
		String[] header = FileUtil.readHeader(resultFile);
		int srcInd = ArrayUtil.indexOf(header, "Source");
		int tgtInd = ArrayUtil.indexOf(header, "Target");
		int relInd = ArrayUtil.indexOf(header, "Relation");
		int srcIDInd = ArrayUtil.indexOf(header, "Source data ID");
		int tgtIDInd = ArrayUtil.indexOf(header, "Target data ID");
		int srcValInd = ArrayUtil.indexOf(header, "Source change");
		int tgtValInd = ArrayUtil.indexOf(header, "Target change");
		int sitesInd = ArrayUtil.indexOf(header, "Sites");

		Set<String> subsetIDs = new HashSet<>();
		Set<String> genes = new HashSet<>();
		Set<String> activatedOncogeneSymbols = new HashSet<>();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + outName + ".sif"));
		Files.lines(Paths.get(resultFile)).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String geneA = t[srcInd];
			String geneB = t[tgtInd];
			String idA = t[srcIDInd];
			String valA = t[srcValInd];
			String idB = t[tgtIDInd];
			String valB = t[tgtValInd];
			boolean upA = !valA.startsWith("-");
			boolean upB = !valB.startsWith("-");

			if (isOncogenic(geneA, idA, upA) && isActivated(geneA, idA, upA))
			{
				subsetIDs.add(idA);
				genes.add(geneA);
				activatedOncogeneSymbols.add(geneA);
			}
			if (isOncogenic(geneB, idB, upB) && isActivated(geneB, idB, upB))
			{
				subsetIDs.add(idB);
				genes.add(geneB);
				activatedOncogeneSymbols.add(geneB);
			}

			if (select(t, geneexp, srcInd, relInd, tgtInd, srcIDInd, srcValInd, tgtIDInd, tgtValInd))
			{
				FileUtil.writeln(ArrayUtil.merge("\t", t[srcInd], t[relInd], t[tgtInd], "", t[sitesInd]), writer);
				subsetIDs.add(idA);
				subsetIDs.add(idB);
				genes.add(geneA);
				genes.add(geneB);
			}
		});

		activatedOncogeneSymbols.stream().sorted().forEach(g -> FileUtil.writeln(g, writer));

		writer.close();

		writeSubsetFormat(dir + "/causative.format", dir + outName + ".format", subsetIDs, genes);
	}

	public static void writeSubsetFormat(String inFile, String outFile, Set<String> ids, Set<String> genes) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("node\tall-nodes\tcolor\t255 255 255\nnode\tall-nodes\tbordercolor\t50 50 50");

		getBackgroundColors(inFile, ids).forEach(l -> FileUtil.lnwrite(l, writer));

		Files.lines(Paths.get(inFile)).forEach(l ->
		{
			String[] t = l.split("\t");

			if (t[2].equals("rppasite") && ids.contains(t[3].split("\\|")[0]) || (genes.contains(t[1]) && (t[3].contains("|!|") || t[3].contains("|i|"))))
			{
				FileUtil.lnwrite(l, writer);
			}
			if (!t[2].equals("rppasite") && ids.contains(t[1]))
			{
				FileUtil.lnwrite(l, writer);
			}
		});

		for (String gene : genes)
		{
			String color = OncoKB.get().isOncogeneOnly(gene) ? "255 50 50" : OncoKB.get().isTumorSuppressorOnly(gene) ? "50 50 255" : OncoKB.get().isCancerGene(gene) ? "180 180 50" : null;

			if (color != null) FileUtil.lnwrite("node\t" + gene + "\tbordercolor\t" + color, writer);
		}

//		if (goi != null) goi.forEach(g -> FileUtil.writeln("node\t" + g + "\tborderwidth\t2", writer));

		writer.close();
	}

	private static Set<String> getBackgroundColors(String inFile, Set<String> ids)
	{
		Set<String> bgColorLines = new HashSet<>();

		Map<String, String> idToSym = FileUtil.linesTabbed(inFile).filter(t -> t[2].equals("tooltip")).collect(Collectors.toMap(t -> t[3].substring(0, t[3].indexOf(",")), t -> t[1]));
		Set<String> symsToCare = ids.stream().filter(idToSym::containsKey).map(idToSym::get).collect(Collectors.toSet());
		FileUtil.lines(inFile).forEach(l ->
		{
			String[] t = l.split("\t");
			if (t.length > 3 && symsToCare.contains(t[1]) && (t[2].equals("color") || t[2].equals("tooltip")))
			{
				bgColorLines.add(l);
			}
		});
		return bgColorLines;
	}

	public static void writeSubsetFormatCorrelation(String inFile, String outFile, Set<String> ids, Set<String> genes) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("node\tall-nodes\tcolor\t255 255 255\nnode\tall-nodes\tbordercolor\t50 50 50");

		Files.lines(Paths.get(inFile)).forEach(l ->
		{
			String[] t = l.split("\t");

			if (t[2].equals("rppasite") && ids.contains(t[3].split("\\|")[0]))
			{
				FileUtil.lnwrite(l, writer);
			}
			if (!t[2].equals("rppasite") && genes.contains(t[1]))
			{
				FileUtil.lnwrite(l, writer);
			}
		});

		for (String gene : genes)
		{
			String color = OncoKB.get().isOncogeneOnly(gene) ? "255 50 50" : OncoKB.get().isTumorSuppressorOnly(gene) ? "50 50 255" : OncoKB.get().isCancerGene(gene) ? "180 180 50" : null;

			if (color != null) FileUtil.lnwrite("node\t" + gene + "\tbordercolor\t" + color, writer);
		}

		writer.close();
	}

	private static boolean select(String[] t, Boolean geneexp, int srcInd, int relInd, int tgtInd, int srcIDInd,
								  int srcValInd, int tgtIDInd, int tgtValInd)
	{
		String geneA = t[srcInd];
		String rel = t[relInd];
		String geneB = t[tgtInd];
		String idA = t[srcIDInd];
		String valA = t[srcValInd];
		String idB = t[tgtIDInd];
		String valB = t[tgtValInd];
		boolean upA = !valA.startsWith("-");
		boolean upB = !valB.startsWith("-");

//		idA = fixTheID(idA, geneA);
//		idB = fixTheID(idB, geneB);

		return isOncogenic(geneB, idB, upB) && // target change is oncogenic
			!isOncogenic(geneA, idA, !upA) && // reverse change of A is not oncogenic
				(geneexp == null || (geneexp == rel.endsWith("expression"))) && // relation is of interest
					isActivated(geneA, idA, upA); // upstream is activated

	}

	private static boolean isActivated(String gene, String id, boolean up)
	{
		if (id.equals(gene) || id.endsWith("-by-network-sig") || id.equals(gene + "-rna") || id.endsWith("_G"))
		{
			return up;
		}
		else
		{
			for (String site : id.substring(gene.length() + 1).split("[-_]"))
			{
				if (site.length() > 1)
				{
					Integer eff = IdentifyOncogenicProteomicChanges.getSiteEffect(gene, site);
					if (eff != null && eff != 0)
					{
						return eff > 0 == up;
					}
				}
			}
		}
		return false;
	}

	/**
	 * For fixing the CPTAC-PanCAn dataset IDs for phosphoproteins.
	 */
	private static String fixTheID(String id, String gene)
	{
		if (!id.startsWith(gene) && id.contains("_"))
		{
			String siteStr = id.split("_")[2];
			String[] sites = siteStr.split("s|t|y");
			id = gene;
			for (String site : sites)
			{
				id += "-" + site;
			}
		}
		return id;
	}

	private static boolean isOncogenic(String gene, String id, boolean up)
	{
		if (id.equals(gene) || id.endsWith("-by-network-sig") || id.equals(gene + "-rna") || id.endsWith("_G") ||  id.endsWith("_R") || id.endsWith("_C"))
		{
			return IdentifyOncogenicProteomicChanges.isChangeOncogenic(gene, null, null, up);
//			return IdentifyOncogenicProteomicChanges.isChangeOncogenicForLung(gene, null, null, up);
		}
		else
		{
			Feature mod = id.endsWith("-P") || id.contains("-P-") || id.endsWith("_P") ? Feature.PHOSPHORYLATION : id.endsWith("-A") || id.contains("-A-") || id.endsWith("_A")? Feature.ACETYLATION : id.endsWith("_C") ? Feature.METABOLITE : null;
			if (mod == null) throw new RuntimeException("Modification is unknown. ID: " + id);

			for (String site : id.substring(gene.length() + 1).split("[-_]"))
			{
				if (site.length() > 1 && IdentifyOncogenicProteomicChanges.isChangeOncogenic(gene, site, mod, up))
//				if (site.length() > 1 && IdentifyOncogenicProteomicChanges.isChangeOncogenicForLung(gene, site, mod, up))
				{
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * Outdated method. Result file structure has changed.
	 */
	private static void writeOncogenicRelationsForCorrelation(String dir) throws IOException
	{
		String outName = "/oncogenic-changes";

		Set<String> subsetIDs = new HashSet<>();
		Set<String> genes = new HashSet<>();

		SiteEffectCollective sec = new SiteEffectCollective();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + outName + ".sif"));
		Files.lines(Paths.get(dir + "/results.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String geneA = t[0];
			String geneB = t[2];
			String idA = t[4];
			String idB = t[5];
			String sites = t[3];

			if (OncoKB.get().isCancerGene(geneB))
			{
				SignedType edgeType = SignedType.typeOf(t[1]);
				boolean select = !edgeType.isSiteSpecific();

				if (!select)
				{
					for (String site : sites.split(";"))
					{
						Feature mod = site.startsWith("K") ? Feature.ACETYLATION : Feature.PHOSPHORYLATION;
						Integer effect = sec.getEffect(geneB, site, mod);
						if (effect != null && effect != 0)
						{
							select = true;
							break;
						}
					}
				}

				if (select)
				{
					FileUtil.writeln(ArrayUtil.merge("\t", t[0], t[1], t[2], "", t[3]), writer);
					subsetIDs.add(idA);
					subsetIDs.add(idB);
					genes.add(geneA);
					genes.add(geneB);
				}
			}
		});

		writer.close();

		writeSubsetFormatCorrelation(dir + "/causative.format", dir + outName + ".format", subsetIDs, genes);

	}
}
