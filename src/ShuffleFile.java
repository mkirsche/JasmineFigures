import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;

/*
 * Shuffles the lines in a file
 */
public class ShuffleFile {
	
	static String fn = "", ofn = "";
	static char header_char = (char)0;
	static int seed = -1;
	static void usage()
	{
		System.out.println();
		System.out.println("Usage: java -cp src ShuffleFile [args]");
		System.out.println("  Example: java -cp src ShuffleFile in_file=in.txt out_file=out.txt");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  in_file          (String) - the file to be shuffled");
		System.out.println("  out_file         (String) - where to write the shuffled file to");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  header_character (char)   - the character indicating header lines");
		System.out.println("  seed             (int)    - the random seed to use when shuffling the lines");
		System.out.println();
	}
	static void parseArgs(String[] args)
	{
		for(String arg : args)
		{
			int equalsIdx = arg.indexOf('=');
			if(equalsIdx == -1)
			{
			}
			else
			{
				String key = arg.substring(0, equalsIdx);
				String val = arg.substring(1 + equalsIdx);
				if(key.equalsIgnoreCase("in_file"))
				{
					fn = val;
				}
				else if(key.equalsIgnoreCase("out_file"))
				{
					ofn = val;
				}
				else if(key.equalsIgnoreCase("header_character"))
				{
					 header_char = val.charAt(0);
				}
				else if(key.equalsIgnoreCase("header_character"))
				{
					 seed = Integer.parseInt(val);
				}
			}
		}
		if(fn.length() == 0 || ofn.length() == 0)
		{
			usage();
			System.exit(0);
		}
	}
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		PrintWriter out = new PrintWriter(new File(ofn));
		
		ArrayList<String> lines = new ArrayList<String>();
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(header_char > 0 && line.charAt(0) == header_char)
			{
				out.println(line);
			}
			else
			{
				lines.add(line);
			}
		}
		
		if(seed == -1)
		{
			Collections.shuffle(lines);
		}
		else
		{
			Collections.shuffle(lines, new Random(seed));
		}
		
		for(String line : lines)
		{
			out.println(line);
		}
		
		input.close();
		out.close();
	}
	
}
