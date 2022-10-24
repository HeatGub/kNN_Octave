%Supervised ML algorithm - k Nearest Neighbours (kNN) with merge sorting for Octave.
%The data is divided equally for training and testing sets.
%Verification of effectiveness is in the form of a plot of draws and successfully assigned classes as a function of k (number of nearest neighbours used to classify the record).
%my dataset was CTG.xls - https://archive.ics.uci.edu/ml/datasets/Cardiotocography - i chose the NSP to be the objective variable (last column: Normal=1; Suspect=2; Pathologic=3). 
%Remember to clean the data from NaNs, empty cells and headings - like in the CTGdata100.xlsx

clear all % remove old variables

%LOADING DATA FROM FILE
            tic
dataraw = xlsread("CTGdata100.xlsx"); % pkg load io (In command window. Also select the directiory where you have the input file)
data = dataraw(:,1:end-2); %remove two last columns in the dataset - objective variable(s) - to be found by the algorithm
lengthR = length(data(:,1)); %rows
lengthC = length(data(1,:)); %columns

% RANDOMIZING INPUT DATA - 50/50 training to testing 
permutt = randperm(lengthR);
dtki = permutt(1:round(lengthR/2)); % KnownDataIndex = Testing data
dtui= permutt(round(lengthR/2)+1:lengthR); % UnknownDataIndex = Training data

%WAITBAR FOR DISTANCES
wtbr = waitbar (0,sprintf("Calculating (distances)^2: %d  ", 0));

% CALCULATING THE SQUARES OF DIFFERENCES (always a positive value)
differsq = zeros(lengthR,lengthR);
for x=1:lengthR
  for i=1:lengthR
    if i==x continue
    elseif any (dtki == i) continue % skip dtki
    elseif any (dtui == x) continue % skip dtui
    endif
      for j = 1:lengthC
          differsq (i,x) = differsq (i,x) + (data(i,j) - data(x,j))^2;
      endfor
  endfor
  waitbar (x/lengthR,wtbr,sprintf("Calculating (distances)^2: %d %% ", round(100*x/lengthR)));
endfor
delete(wtbr)

% CALCULATING THE DISTANCES (square root of the differsq)
Dists(lengthR,lengthR) = 0;
for i=1:lengthR
  if any (dtki == i) continue % skip dtki
  endif
     for j=1:lengthR
      Dists(i,j) = sqrt(differsq(i,j));
  endfor
endfor
lengthRD = length(Dists(:,1));
lengthCD = length(Dists(1,:));
            toc
            tic
           
%INDEXES MATRIX FOR SORTING
MtrxIndxs=zeros(lengthR,lengthR);
for j=1:lengthR
  MtrxIndxs(j,:) = 1:lengthR;
endfor

%WAITBAR FOR SORTING
skip = 0;
wb = waitbar (0,sprintf("sorting: %d / %d ", j-skip, lengthRD));

%MERGE SORT - BEGINNING
function A = mergeswap(A, i, j)
    tmp = A(i);
    A(i) = A(j);
    A(j) = tmp;
end
for j=1:lengthRD
  if any (dtki == j) skip++; continue % skip dtki
  endif
  waitbar (j/lengthRD,wb,sprintf("sorting: %d / %d ", j-skip, round(lengthRD/2)));
  function [D,I] = recmerge(Dists, MtrxIndxs)
      if length(Dists) == 1
          D = Dists;
          I = MtrxIndxs;
      elseif length(Dists) == 2
          if Dists(2) < Dists(1)
              D = mergeswap(Dists, 1, 2);
              I = mergeswap(MtrxIndxs, 1, 2);
          else
              D = Dists;
              I = MtrxIndxs;
          endif
      else
          mid = floor(length(Dists) / 2);
          fin = length(Dists);
          [A,AA]= recmerge(Dists(1:mid),MtrxIndxs(1:mid));
          [B,BB] = recmerge(Dists(mid+1:end),MtrxIndxs(mid+1:end));
          C = [A B];
          CC = [AA BB];
          a = 1;
          b = mid+1;
          for i = 1:fin
              if a < mid+1 && (b >= fin+1 || C(a) <= C(b))
                  D(i) = C(a);
                  I(i) = CC(a);
                  a = a + 1;
              else
                  D(i) = C(b);
                  I(i) = CC(b);
                  b = b + 1;
              endif
          endfor
      endif
  endfunction
  [Dists(j,:),MtrxIndxs(j,:)] = recmerge(Dists(j,:),MtrxIndxs(j,:));
endfor
%MERGE SORT - OVER
delete(wb)
            toc
            tic

%VOTING OF THE k NEIGHBOURS (k is variable to make a plot at the end)
Vote (lengthRD,3) = 0;
for k = 1:50      %k of nearest neighbours
  klist(k)=k;
  clear Vote
  Vote (lengthRD,3) = 0;
  for j=1:lengthRD
    if any (dtki == j) continue % skip dtki
    endif
      for i=round(lengthCD/2):round(lengthCD/2)-1+k %(k = +1:lim)
        if dataraw(MtrxIndxs(j,i),end) == 1 
          Vote(j,1)++; %vote for option 1
        elseif dataraw(MtrxIndxs(j,i),end) == 2
          Vote(j,2)++; %vote for option 2
        elseif dataraw(MtrxIndxs(j,i),end) == 3
          Vote(j,3)++; %vote for option 3
        endif
      endfor
  endfor

  %ANY DRAWS IN THE VOTING?
  draw(k)=0;
  clear draw
  for i=1:lengthR
  maxval(i) = max(Vote(i,:));
  endfor
  for i=1:lengthR
  sumd(i) = sum (Vote(i,:) == maxval(i));
    if any (dtki == i) 
    draw(i)= 0;
    continue % skip dtki 
    elseif  sumd (i) >= 2;
    draw(i)=1;
    endif
  endfor
  draws(k)=length(nonzeros(draw));
  drawspercent(k)=100*(draws(k)/length(dtui));
    
  %ASSIGN THE CLASS BY THE VOTINGS
  Class (lengthR,1) = 0;
  for i=1:lengthR
    if any (dtki == i) continue % skip dtki
    endif
    [maxi,col] = max(Vote(i,:));
    if col == 1
    Class(i)=col;
    elseif col == 2
    Class(i)=col;
    elseif col == 3
    Class(i)=col;
    endif
  endfor
        
  %VERIFICATION
  err = 0;
  for i=1:lengthR
    if any (dtki == i) continue % skip dtki
    endif
    if Class(i) != dataraw(i,end)
    err++;
    endif
  endfor

%EFFECTIVENESS - PLOT
Efctv = 100*((length(dtui)-err)/length(dtui));
Efctvns(k)= Efctv;
endfor
kNNplot = plot(klist,Efctvns,".","markersize", 15, drawspercent,"*","markersize", 5,"color",'r');
grid on, grid minor
xlabel ("k (number of neighbours)");
ylabel ("[%]");
legend('correctly assigned classes','draws in voting',"location",'east');
            toc