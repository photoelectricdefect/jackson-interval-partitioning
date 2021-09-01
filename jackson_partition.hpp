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
     * @param  ptr: Pointer for passing extra data
     * @retval Returns indexes of partitions in data vector
     */
    std::vector<int> partition(std::vector<T> data,fitness_function<T> fn_fit,void* ptr=nullptr) {
      int N=data.size();
      std::vector<int> lastchange(N);
      std::vector<double> opt(N+1);
      opt[0]=0;

      for (size_t i=0;i<N;i++)
      {
        int idx_change=0;
        double opt_cost=-DBL_MAX;

        for (size_t j=0;j<=i;j++)
        {
          double cost=opt[j]+fn_fit(data,j,i,ptr);

          if(cost>opt_cost) {
            opt_cost=cost;
            idx_change=j;
          }
        }

        opt[i+1]=opt_cost;
        lastchange[i]=idx_change;
      }

      std::vector<int> partitions;
      int p=N;

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