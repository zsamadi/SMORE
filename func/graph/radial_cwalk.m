function node_list=radial_cwalk(gStruct, N,wlen, isCheckRadial)

if nargin<4
    isCheckRadial=true;
end
% random walk where each point is progressively farther from the starting point

%     """This function returns sampled paths from a graph.
%     It makes sure points get increasing farther from
%     the starting point.
%     Input:
%     G: the graph
%     N: number of iterations
%     P: probability of ending the walk at each step (%)"""
    XYCoords=gStruct.coordinates;

    node_list = cell(N, 1);
    m=1;

    while m<=N

                %pick a random starting node
                
        start = randi(length(gStruct.neighs));

        %add the index and cell type of the starting node to lists
        temp_node=zeros(1,wlen);
        temp_node(1) = start;
        temp_dist = 0;

        %check if the node has any neighbors
        neighbors = gStruct.neighs{start};
        iWlk=2;


        if ~isempty(neighbors) 

            while iWlk<=wlen && rand>0.2
                neighbors =  gStruct.neighs{start};
                w = gStruct.weights{start};
                probs=w/sum(w);
                next_nodei =randsample(length(probs),1,true,probs) ;
                next_node=neighbors(next_nodei);
                if isCheckRadial
                    next_dist = norm(XYCoords(next_node, :)-XYCoords(temp_node(1), :));
                    isAccept=next_dist > temp_dist;
                else
                    next_dist=0;
                    isAccept=~any(temp_node==next_node);
                end
                
                if isAccept
%                 if next_dist > 0
                    temp_node(iWlk)=next_node;
                    start = next_node;
                    temp_dist = next_dist;
                    iWlk=iWlk+1;
                end
            end
        end
      if iWlk==wlen+1
          if isCheckRadial
             rvDist=XYCoords(temp_node(1:end-1), :)-XYCoords(temp_node(end), :);
             rvDist=sum(rvDist.^2, 2);
             isAccept=all(rvDist(1:end-1)>=rvDist(2:end));
          else

             isAccept=true;
          end

         if isAccept
               node_list{m}=temp_node;
               m=m+1;
         end
      end
      
    end


