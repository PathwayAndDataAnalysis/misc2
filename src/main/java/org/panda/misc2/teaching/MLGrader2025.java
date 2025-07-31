package org.panda.misc2.teaching;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.TermCounter;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 *
 */
public class MLGrader2025
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/ozgunbabur/Documents/Teaching/MachineLearning/Spring2025/";
		Map<String, Double> grades = getOverallGrades(
			dir + "2025-05-23T1117_Grades-CS_438_P_1_01__Applied_Machine_Learning_Spring_2025.csv",
			dir + "CS_438_638_Spring_2025_grades.csv");

		TermCounter tc = new TermCounter();
		grades.keySet().stream().sorted(Comparator.comparing(grades::get))
			.peek(name -> tc.addTerm(getLetterGrade(grades.get(name))))
			.forEach(name -> System.out.println(getLetterGrade(grades.get(name)) + "\t" + grades.get(name) + "\t" + name));
		System.out.println();
		tc.print();
	}

	private static Map<String, Double> getOverallGrades(String canvasFile, String gradescopeFile) throws IOException
	{
		Map<String, Double> gradesMap = new HashMap<>();

		Map<String, double[]> epGrades = readExamsAndProject(canvasFile);
		Map<String, Double> quizGrades = readGroupGradeAvgCanvas(canvasFile, "Quiz ");
		Map<String, Double> paGrades = readGroupGradeAvgGradescope(gradescopeFile, "Programming Assignment ");
		Map<String, Double> hwGrades = readGroupGradeAvgGradescope(gradescopeFile, "Homework ");

		List<String> names = epGrades.keySet().stream().sorted().collect(Collectors.toList());

		for (String name : names)
		{
			if (!paGrades.containsKey(name))
				System.err.println("Missing on gradescope: " + name);
		}

		for (String name : names)
		{
			double exam1 = epGrades.get(name)[0];
			double exam2 = epGrades.get(name)[1];

//			if (exam1 == 0) exam1 = exam2;
//			if (exam2 == 0) exam2 = exam1;

			double project = epGrades.get(name)[2];
			double quiz = quizGrades.get(name);
			double pa = paGrades.get(name);
			double hw = hwGrades.get(name);

			double overall = (0.2 * exam1) + (0.2 * exam2) + (0.05 * quiz) + (0.05 * hw) + (0.3 * project) + (0.2 * pa);
			gradesMap.put(name, overall);
		}

		shiftCurveBy(gradesMap, 100 - 93.97239583333334 + (25.0/29) + 0.91);

		return gradesMap;
	}

	private static Map<String, Double> readGroupGradeAvgCanvas(String file, String prefix)
	{
		Map<String, Double> avgMap = new HashMap<>();
		Map<String, List<Double>> grdMap = new HashMap<>();
		String[] header = FileUtil.readFirstLine(file).split(",");

		FileUtil.lines(file).skip(3).map(l -> l.split(",")).forEach(t ->
		{
			String name = t[0].substring(1) + ", " + t[1].substring(1, t[1].length()-1);
			List<Double> list = new ArrayList<>();
			grdMap.put(name, list);

			for (int i = 5; i < t.length; i++)
			{
				if (t[i].isEmpty()) t[i] = "0";
				if (header[i].startsWith(prefix) && !t[i].equals("NA"))
				{
					list.add(Double.parseDouble(t[i]));
				}
			}
		});

		grdMap.forEach((name, grds) -> avgMap.put(name, CollectionUtil.averageInListOfDouble(grds)));
		return avgMap;
	}

	private static Map<String, Double> readGroupGradeAvgGradescope(String gradescopeFile, String prefix)
	{
		Map<String, Double> avgMap = new HashMap<>();
		Map<String, List<Double>> grdMap = new HashMap<>();
		String[] header = FileUtil.lines(gradescopeFile).findFirst().get().split(",");

		FileUtil.lines(gradescopeFile).skip(1).map(l -> l.split(",")).forEach(t ->
		{
			String name = t[1] + ", " + t[0];
			List<Double> list = new ArrayList<>();
			grdMap.put(name, list);

			for (int i = 5; i < t.length; i++)
			{
				if (header[i].startsWith(prefix) && !header[i].contains(" - ") && !t[i].equals("NA"))
				{
					if (t[i].isEmpty()) t[i] = "0";
					list.add(Double.parseDouble(t[i]));
				}
			}
		});

		grdMap.forEach((name, grds) -> avgMap.put(name, CollectionUtil.averageInListOfDouble(grds)));
		return avgMap;
	}


	private static Map<String, double[]> readExamsAndProject(String file)
	{
		Map<String, double[]> grdMap = new HashMap<>();
		String[] header = FileUtil.readFirstLine(file).split(",");

		FileUtil.lines(file).skip(3).map(l -> l.split(",")).forEach(t ->
		{
			String name = t[0].substring(1) + ", " + t[1].substring(1, t[1].length()-1);
			double[] grades = new double[3];
			grdMap.put(name, grades);

			for (int i = 5; i < t.length; i++)
			{
				if (t[i].isEmpty()) t[i] = "0";

				if (header[i].startsWith("Exam I "))
				{
					grades[0] = Double.parseDouble(t[i]);
				}
				else if (header[i].startsWith("Exam II "))
				{
					grades[1] = Double.parseDouble(t[i]);
				}
				else if (header[i].startsWith("The Project "))
				{
					grades[2] = t[i].startsWith("N") ? 0 : Double.parseDouble(t[i]);
				}
			}
		});

		return grdMap;
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
//		if (v < 12) return "FAI";
//		else if (v < 21) return "CAU";
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
