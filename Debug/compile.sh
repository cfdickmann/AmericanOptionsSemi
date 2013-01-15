g++ ../quellen/*.c* -lglpk -lpthread -o AmericanOptionsSemi 


#liste=$(ls ../quellen/*.c*)

#for i in $liste
#do
#g++ -c ../quellen/$i -lpthread -lglpk & disown
#done

#g++ *.o -lpthread -lglpk -o AmericanOptionsSemi
