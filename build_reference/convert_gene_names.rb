gene_file = "./build_reference/library.gene_description.txt"
xm_file = "./build_reference/library.loc_to_XM.txt"
lg_file = "./build_reference/library.loc_to_linkage.txt"
node_file = ARGV[0]
output_file = ARGV[1] + ".reference.txt"

products = Hash.new
file = File.open(gene_file, "r")

colored = Array.new

#Reads information from the gene_description library and relates each loc ID to its description
while(!file.eof)
  line = file.readline
  line = line.chomp
  
  #stores the information on a specific line if it is correct, otherwise the line is printed 
  if((/\| (.*) \((.*)\),/ =~ line) != nil)
    products[$2] = $1
  else
    puts line
  end
  
end


file.close
xm_file_r = File.open(xm_file, "r")

#Reads the library file to create a mapping between transcript accession numbers and locIDs
while(!xm_file_r.eof)
  
  line = xm_file_r.readline
  line = line.chomp
  arr = line.split("\t\t")
  products[arr[0]]  = arr[1] + "\t" + products[arr[0]]
  
end
xm_file_r.close


lg_file_r = File.open(lg_file, "r")
lg_file_r.readline
#Reads the library file to create a mapping between linkage groups and locIDs
while(!lg_file_r.eof)

  line = lg_file_r.readline
  line = line.chomp
  arr = line.split(" ")
  arr[0].split(",").each{|x|
    if(products.has_key? x)
      products[x]  = arr[1] + "\t" + products[x]
    else
     # puts x
    end
  }
end

#loop done to specify that a certain gene has no LG associated with it, this prevents jagged arrays in the future
products.each_pair{|x,y|
  
  if(y.split("\t").length < 3)
    products[x] = "NO LG\t" + y;
  end
    
}

lg_file_r.close

node_file_r = File.open(node_file, "r")
node_file_r.readline

#Reads in the file provided by the user and creates a hash representation of the data
while(!node_file_r.eof)

  line = node_file_r.readline
  line = line.chomp
  arr = line.split("\t")

  x = arr[1].delete '"'
  y = arr[3].delete '"'

  if(products.has_key? x)
    products[x]  = products[x] + "\t" + y 
    colored.push x
  else
    #puts x
  end
end



out = File.open(output_file, "w")
total = 4
out.puts("LOC ID\tLINKAGE GROUP\tTRANSCRIPT ACCESSION IDENTIFIERS\tGENE NAME\tMODULE COLOR")
#outputs information stored in the hash to a user specified file that will now hold a reference for this file
products.each_pair{|key,value|

  if(colored.include? key)
    out.puts key + "\t" + value

  end
}
out.close
