start check_test(test_name, computed, correct, tol=1E-12);
   failed = 0;
   maxDiff = .;
   if nrow(computed) ^= nrow(correct) | ncol(computed) ^= ncol(correct) then
      failed = 1;
   else if any(computed = .) then
      failed = 1;
   else do;
      maxDiff = max(abs(computed - correct));
      failed = (maxDiff > tol);
   end;
   if failed then do;
      msg = cat("--- ", test_name, " --- FAILS ---");
      print msg[L=""], maxDiff computed correct;
   end;
finish;