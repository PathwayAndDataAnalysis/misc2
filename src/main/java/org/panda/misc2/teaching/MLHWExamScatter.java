package org.panda.misc2.teaching;

import org.panda.utility.FileUtil;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MLHWExamScatter
{
	public static void main(String[] args) throws IOException
	{
		scatterPlot();
	}

	public static void scatterPlot() throws IOException
	{
		String dir = "/home/ozgunbabur/Documents/Teaching/MachineLearning/Spring2025/";
		Map<String, Double> hwGrades = readHWGrades(dir + "gradescope.csv");
		Map<String, Double> exGrades = readExamGrades(dir + "canvas.csv");

		for (String name : hwGrades.keySet())
		{
			double hw = hwGrades.get(name);
			if (!exGrades.containsKey(name))
			{
				System.out.println("name not on canvas = " + name);
			}
			double ex = exGrades.get(name);
			System.out.println(name + "\t" + hw + "\t" + ex);
		}
	}

	public static Map<String, Double> readHWGrades(String gradescopeFile) throws IOException
	{
		String[] header = Files.lines(Paths.get(gradescopeFile)).findFirst().get().split(",");
		List<Integer> inds = new ArrayList<>();
		for (int i = 0; i < header.length; i++)
		{
			if (header[i].startsWith("Homework ") && !header[i].contains("-")) inds.add(i);
		}

		System.out.println("inds.size() = " + inds.size());
		Map<String, Double> map = new HashMap<>();

		Files.lines(Paths.get(gradescopeFile)).skip(1).map(l -> l.split(",")).forEach(t ->
		{
			String name = t[0];

			int sum = 0;

			for (Integer ind : inds)
			{
				sum += t[ind].isEmpty() ? 0 : Integer.parseInt(t[ind]);
			}

			map.put(name, sum / (double) inds.size());
		});

		return map;
	}

	public static Map<String, Double> readExamGrades(String canvasFile) throws IOException
	{
		Map<String, Double> map = new HashMap<>();
		Files.lines(Paths.get(canvasFile)).skip(2).map(l -> l.split(",")).forEach(t ->
		{
			String name = t[1].replaceAll("\"", "").trim() + " " + t[0].replaceAll("\"", "").trim();
			double grade = t[6].isEmpty() ? 0 : Double.parseDouble(t[6]);
			map.put(name, grade);
		});
		return map;
	}

}
