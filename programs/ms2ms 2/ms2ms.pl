#!/usr/bin/perl -w
use strict;
use warnings;
########################################################
# Purpose of the program:                                                                                                  #
# It makes  MSA(microsatellite analyser) input file from MS output  data                              #
########################################################

# Reading all files from the user input folder and writing all output files with name like "out_inputfile"
#  in a the same folder  

print "Enter path of your input files folder\n";
my $dir = <>;
chomp $dir;
my @files ;
opendir( DIRHANDLE, $dir ) or die "couldn't open  : $!";
while  ( defined ( my $file = readdir(DIRHANDLE)) ) 
    {
        if (  $file =~ /^\w+/ ) 
        	{  
        		push (@files,$file)
       		}
	}
foreach my $file (@files)
	{   
		my $infile = "$dir/$file";
	    my $outfile = "$dir/out_$file";
	    program ($infile,$outfile);
	}


###################
#  The main SUBROUTINE     #
###################


sub program
{
      my ($file,$outfile)=@_;
my %msa;


#print "Enter your MS program output file name \n";
#my $file  = <>;
#chomp $file;
open ( IN, $file );
my @lines_list1 = <IN>;

my $total_num_populations = 0;
my $total_num_indsamps = 0;
my @array_of_subpup_sizes=Read_firstline_MSFile(\$total_num_populations,\$total_num_indsamps,@lines_list1)  ;

my @lines_list=Modify_positionlines_MSFile(\$total_num_populations,@lines_list1)  ;

my $num_subpopu = scalar @array_of_subpup_sizes;

 #my $outfile ="out_$file";
 open ( OUT ,">$outfile" ); 
 Caliculate_allele  (\%msa,$num_subpopu,$total_num_populations,@lines_list);


 my %msa2; 
 my $pupu_line_start= 1;
 my $pupu_line_end= 1;
 for (my $i=0; $i<= $#array_of_subpup_sizes; $i++)
 					{   
 						 $pupu_line_end =($array_of_subpup_sizes[$i]+$pupu_line_start) -1;
 						 for (my $x=$pupu_line_start; $x<=$pupu_line_end; $x++ ) 
 								{   
 										for (my $samples_count=1; $samples_count<=$total_num_indsamps; $samples_count++ )
 											  {
 													my $k= int (($x+1)/2);
  											 		my $l=0;
 											 		if ($x%2==0)
 											 			{$l=2;}
 											 			else {$l=1;}
 											 		my $ii = $i+1 ;
 											 		$msa2{$ii}{$k}{$samples_count}{$l}{score}=$msa{$samples_count }{$x}{score};
 											 }					
 							
 								} 
 						$pupu_line_start =$pupu_line_end +1;
 					}





 print OUT "2\t\t\t";
 for (my $j =1; $j <= $total_num_indsamps;$j++)
 			{
 				print OUT"1\t\t";	
 			}
 print OUT"\n\t\t\t";
 for (my $j =1; $j <= $total_num_indsamps;$j++)
 			{
 			print OUT"0\t\t";	
 			}
 			
 print OUT"\n\t\t"; 
 
 for (my $j =1; $j <= $total_num_indsamps;$j++)
 			{      
 					print OUT"\t$j.1\t$j.2";	
 			}
 print OUT"\n";

 foreach my $i (sort{$a <=> $b} keys %msa2)
 		{    
 			foreach my $k (sort {$a <=> $b}keys %{$msa2{$i}})
 			    { 
 						   print OUT "Pop$i\td\t1";
 						   foreach my $samples_count (sort{$a <=> $b} keys %{$msa2{$i}{$k}})
 						   			{   
 						   				foreach my $j (sort {$a <=> $b}keys %{$msa2{$i}{$k}{$samples_count}})
 						   						{
 						   						      print OUT "\t$msa2{$i}{$k}{$samples_count}{$j}{score}";
 						   						} 
 						   			
 						   			}
 						   	print OUT "\n";
 						
 						}
 			
 		
 			}
 


 return 1;
 close IN;
 close OUT; 
 
 
 }#sub program



 
 ##################
 #    other SUBROUTINES     #
 ##################
 
 
 
 # reads the first line of the ms output file 
sub Read_firstline_MSFile
	{
		my ($total_num_populations,$total_num_indsamps,@list) = @_;
		my @sizes;
		my $num_subpopu;
		my @words = split /\s+/, $list[0];
		
		$$total_num_populations = $words[1];
		$$total_num_indsamps  =  $words[2];
		if ($list[0] =~ /-m/  ||  $list[0] =~ /-I/)
		 	{  for (my $i= 0; $i <= $#words; $i++)
    				{
	    				if ($words[$i] eq "-m" or $words[$i] eq "-I")
	      					{  
	    						my $position =  $i+1;
	    			 			$num_subpopu = $words[$position];
	    						my $total = $position + $num_subpopu;
	        					for (my $j= $position+1   ; $j <= $total; $j++)
	    								{
		           							push (@sizes,$words[$j]);
	            		    			}
	        				}
      			  }   
      	  }
      else { push (@sizes,$words[1]);}
	  
	  return @sizes;
	}
 
 
 

#modify the lines where segregating sites are zero ,, introduces a line "positions:  "
sub Modify_positionlines_MSFile
	{
		
		my ($total_num_populations,@list) = @_;
		my @list1;
		for (my $i= 0; $i <= $#list; $i++)
    			{
    				if ($list[$i] =~ /^segsites: 0/ )
	      				{  	push (@list1,$list[$i] );
		           				push (@list1,"positions: \n" );
		           				for  (my $j=1;$j<=$$total_num_populations;$j++)
		           						{push (@list1,"\n" );}
		           				
		           			}
	        		 else {push (@list1,$list[$i] );}		
      				}   
	     return @list1;
	}






# caliculating the allele count
sub Caliculate_allele
{ 
     my ($msa,$num_subpopu,$total_num_populations,@list) = @_;
     my $total_lines = $#list ;
     my $total_samples_count=0;
     
     my $samples_count = 1;
     my $subpop_count =1;
     for(my $k=0; $k<=$total_lines; $k++) 
    	  {
			 
			if ( $list[$k] =~ /^\/\// ) 
				{     $total_samples_count ++;
						my @segsites = split /\s+/, $list[($k+1)];
					    my $mutation_matrix_columns= $segsites[1];
				        my $mutation_matrix_rows=  $total_num_populations;
				        my $a = $k+3;
				         my $b = $a+($mutation_matrix_rows -1);
				     	 my @array_scores = scores($a,$b,$mutation_matrix_columns,$mutation_matrix_rows,@list);
				     	 for  (my $x=1; $x <= scalar @array_scores; $x++)
					  		{	 $$msa{$total_samples_count}{$x}{score}=  $array_scores[$x-1];
					  		}
				  } 
   		     }  
   
   
   return 1;
  }



           # calculates the scores based on the mutations
			sub  scores
             			 { 
        				 		 my ($sart_line, $end_line,$columns,$mutation_matrix_rows,@list )    = @_;
                    
                                 my @segsites_array = split(' ', $list[$sart_line-2]);
                            
                                 if ($segsites_array[1]==0)
                                      { my @scores_array;
                          
                                         for(my $y=$sart_line; $y <= $end_line ; $y++)
                                            { 
                                            	push (@scores_array,100);
                                             }
                                         return @scores_array;
                                        }
                                else 
                                {
                                	my @random_number = random_number($columns);
                                	my @scores_array;
							     	for(my $y=$sart_line; $y <= $end_line ; $y++) 
								    	{   $list[$y]  =~ s/\s*//g;
								    	 	my @row_array = split('', $list[$y]);
					    	          	    my $score = 100; 
						        		    for(my $d=0; $d <= $#row_array ; $d++)
		    				         		     {    
		    				         		       if ($row_array[$d] eq 1)
		    				         		      		      {   
		    				         		      		      		if ( $random_number[$d] == 1)
		    				         		          					{   $score++;
		    				         		          					 }
		    			         			                         else 
																		 {   $score--;
		    				         		          					 }
		    					      		     				}
		    		    				      		} 
		    		    				      push (@scores_array,$score);
		    		    				      
							      	  } 
							  		
							        return @scores_array;
							        }
							    
							 }  

		
		# generates the random number	
		sub random_number  
						{ 
								my $size =$_[0];
								my @array= (0,1);
								my @array1;
								for(my $i =1 ; $i <= $size ; $i++)
									{ 
   											push (@array1, int rand @array)
    		   						 }
							  return @array1;
							}


exit;
