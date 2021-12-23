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
import java.util.Set;

/**
 * Assumes each proteomic row has one gene.
 */
public class CausalPathShowPathsToCancer
{

	public static void main(String[] args) throws IOException
	{
//		generateRecursivelyComparison("/Users/ozgun/Documents/Analyses/CPTAC-LSCC-3.2/clusters/against-normals", Boolean.TRUE);
//		generateRecursivelyComparison("/Users/ozgun/Documents/Analyses/CPTAC-LSCC-3.2/temp", null);
		generateRecursivelyComparison("/Users/ozgun/Documents/Analyses/CPTAC-PanCan/clusters/against-others-diffexp", null);
	}

	public static void generateRecursivelyComparison(String dir, Boolean phospho) throws IOException
	{
		if (Files.exists(Paths.get(dir + "/results.txt")))
			writeOncogenicRelations(dir, phospho);

		for (File d : new File(dir).listFiles())
		{
			if (d.isDirectory()) generateRecursivelyComparison(d.getPath(), phospho);
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

	private static void writeOncogenicRelations(String dir, Boolean phospho) throws IOException
	{
		String outName = "/oncogenic-changes" + (phospho == null ? "" : phospho ? "-phospho" : "-express");

		Set<String> subsetIDs = new HashSet<>();
		Set<String> genes = new HashSet<>();
		Set<String> activatedOncogeneSymbols = new HashSet<>();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + outName + ".sif"));
		Files.lines(Paths.get(dir + "/results.txt")).skip(1).map(l -> l.split("\t")).forEach(t ->
		{
			String geneA = t[0];
			String geneB = t[2];
			String idA = t[4];
			String valA = t[5];
			String idB = t[7];
			String valB = t[8];
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

			if (select(t, phospho))
			{
				FileUtil.writeln(ArrayUtil.merge("\t", t[0], t[1], t[2], "", t[3]), writer);
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

		Files.lines(Paths.get(inFile)).forEach(l ->
		{
			String[] t = l.split("\t");

			if (t[2].equals("rppasite") && ids.contains(t[3].split("\\|")[0]) || (genes.contains(t[1]) && t[3].contains("|a|")))
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

	private static boolean select(String[] t, Boolean phospho)
	{
		String geneA = t[0];
		String rel = t[1];
		String geneB = t[2];
		String idA = t[4];
		String valA = t[5];
		String idB = t[7];
		String valB = t[8];
		boolean upA = !valA.startsWith("-");
		boolean upB = !valB.startsWith("-");

		idA = fixTheID(idA, geneA);
		idB = fixTheID(idB, geneB);

		return isOncogenic(geneB, idB, upB) && // target change is oncogenic
			!isOncogenic(geneA, idA, !upA) && // reverse change of A is not oncogenic
				(phospho == null || (phospho == rel.contains("phospho"))) && // relation is of interest
					isActivated(geneA, idA, upA); // upstream is activated

	}

	private static boolean isActivated(String gene, String id, boolean up)
	{
		if (id.equals(gene) || id.endsWith("-by-network-sig") || id.equals(gene + "-rna"))
		{
			return up;
		}
		else
		{
			for (String site : id.substring(gene.length() + 1).split("-"))
			{
				Integer eff = IdentifyOncogenicProteomicChanges.getSiteEffect(gene, site);
				if (eff != null && eff != 0)
				{
					return eff > 0 == up;
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
		if (id.equals(gene) || id.endsWith("-by-network-sig") || id.equals(gene + "-rna"))
		{
			if (IdentifyOncogenicProteomicChanges.isChangeOncogenic(gene, null, null, up))
			{
				return true;
			}
		}
		else
		{
			Feature mod = id.endsWith("-P") || id.contains("-P-") ? Feature.PHOSPHORYLATION : id.endsWith("-A") || id.contains("-A-") ? Feature.ACETYLATION : null;
			if (mod == null) throw new RuntimeException("Modification is unknown. ID: " + id);

			for (String site : id.substring(gene.length() + 1).split("-"))
			{
				if (site.length() > 1 && IdentifyOncogenicProteomicChanges.isChangeOncogenic(gene, site, mod, up))
				{
					return true;
				}
			}
		}
		return false;
	}

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
