function [node_list,length_list]=radial_balanced_cwalk(gStruct,N,min_w, max_w)
%     """This function returns sampled paths from a graph.
%     It makes sure points get increasing farther from
%     the starting point and tries to balance the sampling
%     by adjusting the probability by which an edge is chosen
%     based on how many times it has been traversed before.
%     Input:
%     G: the graph
%     N: number of iterations
%     P: probability of ending the walk at each step (%)
%     edge_labels: unique identifiers of the edges
%     min_w: minimum accepted length for each path"""

if nargin<5
    max_w=min_w;
end


gStruct.neighsCovered=cell(length(gStruct.neighs), 1);

for ig=1:length(gStruct.neighs)
     gStruct.neighsCovered{ig}=ones(length(gStruct.neighs{ig}), 1);
end

        
    XYCoords=gStruct.coordinates;
    node_list = cell(N, 1);
%     edge_list =  cell(N, 1);
    length_list=zeros(N,1);

    m=1;
    while(m<=N) 
        %pick a random starting node
        start = randi(length(gStruct.neighs));

        %add the index and cell type of the starting node to lists
        temp_node=zeros( 1, max_w);
%         temp_edge=zeros(ceil(min_w+(1/P)^2), 2);
        temp_node(1) = start;
        temp_dist = 0;

        %check if the node has any neighbors
        neighbor = gStruct.neighs{start};
        iWlk=2;
        
        if ~isempty(neighbor) 

            while iWlk<=max_w && rand>0.2
                outgoing_edge_index =  gStruct.neighs{start};
                outgoing_edge_counts = gStruct.neighsCovered{start};
                w = 1./outgoing_edge_counts.^4; %4 is chosed based on a couple of trials

                probs=w/sum(w);

                next_nodei =randsample(length(probs),1,true,probs) ;



%                 next_nodei =randsample(length(probs),1,true,probs) ;
                next_node=outgoing_edge_index(next_nodei);
                next_dist = norm(XYCoords(next_node, :)-XYCoords(temp_node(1), :));

                if next_dist > temp_dist
%                 if next_dist > 0

                    temp_node(iWlk)=next_node;
%                     temp_edge(iWlk-1, :)=[start,next_node];
                    temp_dist = next_dist;
                    gStruct.neighsCovered{start}(next_nodei)=gStruct.neighsCovered{start}(next_nodei)+1;
                    start = next_node;
                    iWlk=iWlk+1;
                end
            end
        end
        temp_node=temp_node(1:iWlk-1);
%         temp_edge=temp_edge(1:iWlk-2, :);

        if iWlk>min_w

         rvDist=XYCoords(temp_node(1:end-1), :)-XYCoords(temp_node(end), :);
         rvDist=sum(rvDist.^2, 2);
         if all(rvDist(1:end-1)>rvDist(2:end))
            node_list{m}=temp_node;
    %         edge_list{m}=temp_edge;
            length_list(m)=length(temp_node);
    
            m=m+1;
         end




        end
        
    end
