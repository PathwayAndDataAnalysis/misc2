package org.panda.misc2.teaching;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.TermCounter;
import org.panda.utility.statistics.Summary;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * To fix blackboard file, open it with sublime, cop and paste the text in a new file and save with the same name.
 */
public class DMGraderFall2022
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/ozgunbabur/Documents/Teaching/DiscreteMath/Fall2022/grading/";
		Map<String, Double> grades = getOverallGrades(dir + "CS_220_Fall_2022_grades.csv",
			dir + "gc_B3210-1482_fullgc_2022-12-27-11-38-23.tsv",
			dir + "UMBCS220MATH320BaburFall2022_assignment_report_2022-12-19_1414.csv");

		TermCounter tc = new TermCounter();
		grades.keySet().stream()
			.sorted(Comparator.comparing(grades::get))
//			.sorted()
			.peek(name -> tc.addTerm(getLetterGrade(grades.get(name)))).forEach(name -> System.out.println(getLetterGrade(grades.get(name)) + "\t" + grades.get(name) + "\t" + name));
		System.out.println();
		tc.print();
	}

	private static Map<String, Double> getOverallGrades(String gradescopeFile, String blackboardFile, String zyBookFile) throws IOException
	{
		Map<String, Double> gradesMap = new HashMap<>();

		Map<String, Double> zyBookGrades = readZyBookGrades(zyBookFile, zyBookFile + "-summed.txt");
		Map<String, Double> hwGrades = readHomeworkGradeAvg(gradescopeFile);
		Map<String, double[]> mfGrades = readMidtermAndFinal(blackboardFile);

		List<String> names = hwGrades.keySet().stream().sorted().collect(Collectors.toList());

		for (String name : names)
		{
			if (!zyBookGrades.containsKey(name))
			{
				System.err.println("zyBook misses: " + name);
				zyBookGrades.put(name, 0D);
			}
			if (!hwGrades.containsKey(name))
			{
				System.err.println("Gradescope misses: " + name);
				hwGrades.put(name, 0D);
			}
			if (!mfGrades.containsKey(name))
			{
				System.err.println("Blackboard misses: " + name);
				mfGrades.put(name, new double[]{0, 0});
			}

			double midterm = mfGrades.get(name)[0];
			double finalExam = mfGrades.get(name)[1];
			double hw = hwGrades.get(name);
			double book = zyBookGrades.get(name);

			double overall = (0.3 * midterm) + (0.4 * finalExam) + (0.15 * book) + (0.15 * hw);
//			double overall = (0.5 * book) + (0.5 * hw);
			gradesMap.put(name, overall);
		}

		shiftCurveBy(gradesMap, (100 - 89.29399366578691));// + (55 - 54.1040731629393));

		return gradesMap;
	}

	private static Map<String, Double> readZyBookGrades(String inFile, String outFile) throws IOException
	{
		String[] header = Files.lines(Paths.get(inFile)).findFirst().get().split(",");
		double[] weights = new double[header.length];

		for (int i = 3; i < header.length; i++)
		{
			weights[i] = Double.parseDouble(header[i].substring(header[i].indexOf("(") + 1, header[i].indexOf(")")));
		}

		double totW = Summary.sum(weights);
		for (int i = 0; i < weights.length; i++)
		{
			weights[i] /= totW;
		}

		Map<String, Double> overall = new HashMap<>();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write(header[0] + "\t" + header[1] + "\tGrade");

		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split(",")).forEach(t ->
		{
			String name = t[0] + ", " + t[1];
			FileUtil.lnwrite(t[0] + "\t" + t[1] + "\t", writer);

			double tot = 0;
			for (int i = 3; i < weights.length; i++)
			{
				tot += Double.parseDouble(t[i]) * weights[i];
			}

			overall.put(name, tot);
			FileUtil.write("" + (int) Math.round(tot), writer);
		});

		writer.close();
		return overall;
	}

	private static Map<String, Double> readHomeworkGradeAvg(String file)
	{
		Map<String, Double> avgMap = new HashMap<>();
		Map<String, List<Double>> grdMap = new HashMap<>();
		String[] header = FileUtil.lines(file).findFirst().get().split(",");

		FileUtil.lines(file).skip(1).map(l -> l.split(",")).forEach(t ->
		{
			String[] s = t[0].split(" ");
			String name = (s.length > 1 ? s[1] + ", " : "") + s[0];
			List<Double> list = new ArrayList<>();
			grdMap.put(name, list);

			for (int i = 3; i < t.length; i++)
			{
				if (header[i].startsWith("Homework ") && !header[i].contains(" - ") && !t[i].equals("NA"))
				{
					if (t[i].isEmpty()) t[i] = "0";
					list.add(Double.parseDouble(t[i]));
				}
			}
		});

		grdMap.forEach((name, grds) -> avgMap.put(name, CollectionUtil.averageInListOfDouble(grds)));
		return avgMap;
	}

	private static Map<String, Double> readQuizGradeAvg(String file)
	{
		Map<String, Double> avgMap = new HashMap<>();
		Map<String, List<Double>> grdMap = new HashMap<>();
		String[] header = FileUtil.readHeader(file);

		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			String name = t[0] + ", " + t[1];
			List<Double> list = new ArrayList<>();
			grdMap.put(name, list);

			for (int i = 7; i < t.length; i++)
			{
				if (t[i].isEmpty()) t[i] = "0";
				if (header[i].startsWith("Quiz ") && !t[i].equals("NA"))
				{
					list.add(Double.parseDouble(t[i]));
				}
			}
		});

		grdMap.forEach((name, grds) -> avgMap.put(name, CollectionUtil.averageInListOfDouble(grds)));
		return avgMap;
	}

	private static Map<String, double[]> readMidtermAndFinal(String file)
	{
		Map<String, double[]> grdMap = new HashMap<>();
		String[] header = FileUtil.readHeader(file);

		FileUtil.linesTabbedSkip1(file).forEach(t ->
		{
			String name = t[0] + ", " + t[1];
			double[] grades = new double[2];
			grdMap.put(name, grades);

			for (int i = 3; i < t.length; i++)
			{
				if (t[i].isEmpty()) t[i] = "0";

				if (header[i].startsWith("Midterm"))
				{
					grades[0] = Double.parseDouble(t[i]);
				}
				else if (header[i].startsWith("Final"))
				{
					grades[1] = Double.parseDouble(t[i]);
				}
			}
		});

		return grdMap;
	}

	private static Map<String, Double> readExtraCredit(String file)
	{
		return FileUtil.linesTabbed(file).collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));
	}

	private static void shiftCurveBy(Map<String, Double> grades, double amount)
	{
		for (String name : grades.keySet())
		{
			grades.put(name, grades.get(name) + amount);
		}
	}

//	private static String getLetterGrade(double v)
//	{
//		if (v < 40) return "FAI";
//		else if (v < 70) return "CAU";
//		else return "SAT";
//	}

	private static String getLetterGrade(double v)
	{
		if (v < 40) return "F";
		else if (v < 45) return "D-";
		else if (v < 50) return "D";
		else if (v < 55) return "D+";
		else if (v < 60) return "C-";
		else if (v < 65) return "C";
		else if (v < 70) return "C+";
		else if (v < 75) return "B-";
		else if (v < 80) return "B";
		else if (v < 85) return "B+";
		else if (v < 90) return "A-";
		else return "A";
	}
}
