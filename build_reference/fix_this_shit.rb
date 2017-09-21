#This file serves as a tool for fixing a small issue in one of the library files by fixing the alignment of the file's different components

file = File.open(ARGV[0], "r")

collection = Hash.new

while(!file.eof)

	line = file.readline()
	arr = line.chomp.split("\t")
	

	arr[0].split(",").each{|elem|
		collection[elem] = arr[1];
	}
end

printcollection = Hash.new

collection.each_pair{|key,value|

	printcollection[value] = Array.new if(not printcollection.has_key?(value))
	printcollection[value].push(key);
}


file.close;

out = File.open("output.out","w")
printcollection.each_pair{|key,value|
	
	ret = value.join(",")
	out.write("#{ret}\t#{key}\n")
}
out.close;
