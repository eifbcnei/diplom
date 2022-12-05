function[logger] = log_progress(idx, buffer, total, enabled)
   if (enabled==1)
      percentDone = 100 * idx / total;
      msg = sprintf('Job done: %3.1f//100.0\n', percentDone);
      fprintf([buffer, msg]);
      logger = repmat(sprintf('\b'), 1, length(msg));
   else 
      logger = ''; 
   end
end