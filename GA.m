clear all; clc;

nPop = 100;
nKrom = 6;

[x1,x2,Chromosome] = initializeChromosome();                                %melakukan inialisasi kromosom

for i=1:nPop
    for j=1:nKrom
        fitness(j) = fitnessFunction(x1(j),x2(j));
    end
    probF = probability(fitness);
    cumProb = cumsum(probF);
    
    rwChrom = rouletteWheel(cumProb);
    
    newCrossedChromosome = crossover(Chromosome);
    
    finalChromosome = mutate(newCrossedChromosome);
    
    [newX1,newX2] = DecodeChromosome(finalChromosome);                      %melakukan decode kromosome
    
    [a,b,index] = findBestChromosome(newX1,newX2);                          %melakukan mencari kromosom terbaik
    
    value = objective(a(index),b(index));                                   %memasukkan hasil objektif ke value
    pop{i,1} = cat(2,value,a(index),b(index),finalChromosome(index,:));
    x1 = newX1;
    x2 = newX2;
    Chromosome = finalChromosome;
    clear a b index value newX1 newX2 newCrossedChromosome;
end

pop = cell2mat(pop);
best = min(pop(:,1));
index = find(pop(:,1)==best,1);

bestChromosome = pop(index,4:9);
disp(['Best Chromosome : ',num2str(bestChromosome)]);
disp(['X1 = ',num2str(pop(index,2))]);
disp(['X2 = ',num2str(pop(index,3))]);
disp(['Minimum value : ',num2str(min(pop(:,1)))]);

function [x1,x2,kromosom] = initializeChromosome()

    kromosom = randi([0 1],6,6);
    
    [x1,x2] = DecodeChromosome(kromosom);
end

function [x1,x2] = DecodeChromosome(kromosom)
    bawah1 = -3;                                                            %batas bawah x1
    atas1 = 3;                                                              %batas atas x1
    bawah2 = -2;                                                            %batas bawah x2
    atas2 = 2;                                                              %batas atas x2
    pem = 0;
    for i=1:3
        pem = pem + (2^(-i));
    end
    for j=1:6
        pen1 = 0;
        pen2 = 0;
        nP = 1;
        for k=1:3
            pen1 = pen1 + (kromosom(j,k)*(2^(-nP)));
            nP = nP+1;
        end
        nP = 1;
        for k=4:6
            pen2 = pen2 + (kromosom(j,k)*(2^(-nP)));
            nP = nP+1;
        end
        x1(j) = bawah1 + ((atas1-bawah1)/pem)* pen1;
        x2(j) = bawah2 + ((atas2-bawah2)/pem)* pen2;
    end
end

function o = objective(x1,x2)                                               % rumus fungsi yang dari soal                                            
    pen1 = (4-(2.1*(x1^2))+((x1^4)/3));                                     % pengali yang mengalikan rumus
    pen2 = (pen1*(x1^2))+(x1*x2) + (-4+(4*(x2^2))*(x2^2));                  % melakukan pengalian terhadap pengali
    
    o = pen2;                                                               % return objective ke variabel o
end

function f = fitnessFunction(x1,x2)
    o = objective(x1,x2);

    f = abs(1/(o+0.01));
end

function p = probability(f)
    tot = sum(f);
    
    for i=1:6
        p(i) = f(i) / tot;
    end
end

function newChromosome = rouletteWheel(cumProb)
    roulette = rand(1,6);
    for i=1:6
        for j=1:6
            if roulette(i) < min(cumProb)
                indexChromosome = 0;
                break;
            elseif roulette(i)>cumProb(j) && roulette(i)<cumProb(j+1)
                indexChromosome = j;
                break;
            end
        end
        newChromosome(i) = indexChromosome+1;
    end
end

function newCrossedChromosome = crossover(Chrom)                            % melakukan crossover
    probCross = 0.65;                                                       % probability 
    cross = rand(1,6);
    indexParent=1;
    for i=1:6
        if cross(i)>probCross                                               % melakukan crossovr jika lebih besar dari probabilityCrossOver
            parent(indexParent) = i;
            indexParent = indexParent + 1;
        end
    end
    if exist('parent','var')==0
        newCrossedChromosome = Chrom;
    elseif size(parent,2)==1
        newCrossedChromosome = Chrom;
    else 
        random = randi([2 size(Chrom,1)-1],1);
        for i=1:size(Chrom,1)
            if i ~= size(Chrom,1)
                newCrossedChromosome{i,1} = cat(2,Chrom(i,1:random), Chrom(i+1,random+1:6));
            elseif i == size(Chrom,1)
                newCrossedChromosome{i,1} = cat(2,Chrom(i,1:random), Chrom(1,random+1:6));
            end
        end
        newCrossedChromosome = cell2mat(newCrossedChromosome);
    end
end

function [o,f,index] = findBestChromosome(x1,x2)
    for i=1:6
        o(i) = objective(x1(i),x2(i));
        f(i) = abs(1/(o(i)+0.01));
    end
    bestSoFar = max(f); 
    index=find(f==bestSoFar,1);
end

function chrom = mutate(chrom)                                              % melakukan Mutasi
    mutationProb = 0.065;                                                   % Probability Mutasi
    indexes = rand(size(chrom))<mutationProb;                               % melakukan mutasi jika lebih besar dari Probability Mutasi
    chrom(indexes) = chrom(indexes)*-1+1;
end