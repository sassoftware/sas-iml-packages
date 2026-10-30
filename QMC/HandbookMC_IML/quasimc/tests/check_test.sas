start check_test(test_name, prob, correct, tol=1E-12);
   failed = 0;
   maxDiff = .;
   if nrow(prob) ^= nrow(correct) | ncol(prob) ^= ncol(correct) then
      failed = 1;
   else if any(prob = .) then
      failed = 1;
   else do;
      maxDiff = max(abs(prob - correct));
      failed = (maxDiff > tol);
   end;
   if failed then do;
      msg = cat("--- ", test_name, " --- FAILS ---");
      print msg[L=""], maxDiff prob correct;
   end;
finish;