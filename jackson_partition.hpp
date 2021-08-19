#ifndef JACKSONPT_H_
#define JACKSONPT_H_

#include <vector>
#include <limits>

/**
 * @brief  Fitness function for estimating partition fitness
 * @note   
 * @param  data: Data
 * @param  from: Partition start
 * @param  to: Partition end
 * @param  ptr: Pointer for passing extra data
 * @retval 
 */
template<typename T>
using fitness_function = double(*)(std::vector<T> data,int from,int to,void*ptr);

/**
 * @brief  Class for finding the optimal partitioning of an interval I, given a fitness function which calculates
 * the fitness of each subpartition P
 * @note   
 * @retval None
 */
template<typename T>
class jackson_partition {
  private:
  public:
    /**
     * @brief  Finds optimal partition of data given a partition fitness function in O(N^2)
     * @note   
     * @param  data: Datapoint vector
     * @param  fn_fit: Fitness function
     * @param  minimize: If set to true minimizes fitness function, otherwise it is maximized
     * @param  ptr: Pointer for passing extra data
     * @retval Returns indexes of partitions in data vector
     */
    std::vector<int> partition(std::vector<T> data,fitness_function<T> fn_fit,bool minimize=false,void* ptr=nullptr) {
      int N=data.size();
      std::vector<int> lastchange(N);
      std::vector<double> opt(N);
      opt[0]=0;

      if(minimize) {
        for (size_t i=1;i<=N;i++)
        {
          double opt_cost=DBL_MAX;
          int idx_change=0;

          for (size_t j=1;j<=i;j++)
          {
            double cost=opt[j-1]+fn_fit(data,j,i,ptr);

            if(cost<opt_cost) {
              opt_cost=cost;
              idx_change=j;
            }
          }

          opt[i]=opt_cost;
          lastchange[i]=idx_change;
        }
      }
      else {
        for (size_t i=1;i<=N;i++)
        {
          double opt_cost=-DBL_MAX;
          int idx_change=0;

          for (size_t j=1;j<=i;j++)
          {
            double cost=opt[j-1]+fn_fit(data,j,i,ptr);

            if(cost>opt_cost) {
              opt_cost=cost;
              idx_change=j;
            }
          }

          opt[i]=opt_cost;
          lastchange[i]=idx_change;
        }
      }

      std::vector<int> partitions;
      int p=N+1;

      do {
        p--;
        p=lastchange[p];
        partitions.push_back(p);
      } while(p>0);
    
      std::reverse(partitions.begin(),partitions.end());

      return partitions;
    }
}; 

#endif