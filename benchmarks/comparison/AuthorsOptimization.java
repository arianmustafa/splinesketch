import java.lang.instrument.Instrumentation;
import java.lang.management.ManagementFactory;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.*;

/** Compare a patch to the pinned Java authors' code in the same JVM. */
public class AuthorsOptimization {
    private static Instrumentation instrumentation;
    private static volatile long sink;
    private static long rankChecks, stateChecks, derivativeChecks, mergeChecks;
    private static int sharedFailures;
    public static void premain(String options, Instrumentation inst) { instrumentation = inst; }

    private static final class Sketch {
        final boolean original, mg;
        final Object value;
        Sketch(boolean original, boolean mg, int k) {
            this.original = original; this.mg = mg;
            value = original ? (mg ? new OriginalSplineSketchMG(k) : new OriginalSplineSketch(k))
                             : (mg ? new SplineSketchMG(k) : new SplineSketch(k));
        }
        void update(double x) {
            if (original) { if (mg) ((OriginalSplineSketchMG)value).update(x); else ((OriginalSplineSketch)value).update(x); }
            else { if (mg) ((SplineSketchMG)value).update(x); else ((SplineSketch)value).update(x); }
        }
        void finish() {
            if (original) { if (mg) ((OriginalSplineSketchMG)value).consolidate(); else ((OriginalSplineSketch)value).consolidate(); }
            else { if (mg) ((SplineSketchMG)value).consolidate(); else ((SplineSketch)value).consolidate(); }
        }
        void resize(int k) {
            if (original) { if (mg) ((OriginalSplineSketchMG)value).resize(k); else ((OriginalSplineSketch)value).resize(k); }
            else { if (mg) ((SplineSketchMG)value).resize(k); else ((SplineSketch)value).resize(k); }
        }
        void compress() {
            if (!mg) throw new IllegalStateException("only the MG class exposes final compression");
            if (original) ((OriginalSplineSketchMG)value).compressNonFrequentToBucketsAndResize();
            else ((SplineSketchMG)value).compressNonFrequentToBucketsAndResize();
        }
        List<Integer> query(List<Double> xs) {
            if (original) return mg ? ((OriginalSplineSketchMG)value).query(xs) : ((OriginalSplineSketch)value).query(xs);
            return mg ? ((SplineSketchMG)value).query(xs) : ((SplineSketch)value).query(xs);
        }
    }

    private static final class RandomBits {
        long state;
        RandomBits(long seed) { state = seed; }
        long next() {
            long z = state += 0x9e3779b97f4a7c15L;
            z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
            z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
            return z ^ (z >>> 31);
        }
        long bounded(long n) { return Long.remainderUnsigned(next(), n); }
    }
    private static double[] data(int shape, int seed, int n) {
        RandomBits rng = new RandomBits(0x6a09e667f3bcc909L + seed * 0x100000001b3L + shape);
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            switch (shape) {
                case 0: result[i] = rng.bounded(1000001) / 1000.0 - 500; break;
                case 1:
                    long sum = 0; for (int j = 0; j < 4; j++) sum += rng.bounded(2001);
                    result[i] = (sum - 4000) / 1000.0; break;
                case 2:
                    long choice = rng.bounded(10);
                    result[i] = choice < 7 ? choice % 3 - 1 : (rng.bounded(129) - 64) / 8.0; break;
                case 3:
                    if (rng.bounded(50) == 0) {
                        double sign = (rng.next() & 1) != 0 ? 1 : -1;
                        result[i] = sign * (10000 + rng.bounded(1000));
                    } else result[i] = (rng.bounded(2001) - 1000) / 1000.0;
                    break;
                default: result[i] = i + rng.bounded(3) / 10.0;
            }
        }
        return result;
    }
    private static List<Double> probes(double[] input) {
        double[] sorted = input.clone(); Arrays.sort(sorted);
        List<Double> result = new ArrayList<>();
        result.add(Double.NEGATIVE_INFINITY); result.add(Double.POSITIVE_INFINITY);
        result.add(-0.0); result.add(0.0);
        for (int i = 0; i <= 256; i++) {
            int at = i * (sorted.length - 1) / 256;
            double x = sorted[at];
            result.add(x); result.add(Math.nextDown(x)); result.add(Math.nextUp(x));
            if (at + 1 < sorted.length && x < sorted[at + 1]) result.add(x + (sorted[at + 1] - x) / 2);
            result.add(sorted[0] + (sorted[sorted.length - 1] - sorted[0]) * i / 256);
        }
        return result;
    }
    private static void same(Object a, Object b) {
        if (a == null || b == null) { if (a != b) throw new AssertionError("null state mismatch"); return; }
        if (a.getClass().isArray()) {
            if (!b.getClass().isArray() || Array.getLength(a) != Array.getLength(b)) throw new AssertionError("array size mismatch");
            for (int i = 0; i < Array.getLength(a); i++) same(Array.get(a, i), Array.get(b, i));
        } else if (a instanceof Map<?, ?>) {
            Map<?, ?> x = (Map<?, ?>)a, y = (Map<?, ?>)b;
            if (!x.keySet().equals(y.keySet())) throw new AssertionError("map keys mismatch");
            for (Object key : x.keySet()) same(x.get(key), y.get(key));
        } else if (!a.equals(b)) throw new AssertionError("state mismatch: " + a + " versus " + b);
    }
    private static void state(Sketch a, Sketch b) throws Exception {
        for (Field f : a.value.getClass().getDeclaredFields()) {
            if (Modifier.isStatic(f.getModifiers())) continue;
            Field g = b.value.getClass().getDeclaredField(f.getName());
            f.setAccessible(true); g.setAccessible(true);
            try { same(f.get(a.value), g.get(b.value)); }
            catch (AssertionError error) { throw new AssertionError(f.getName(), error); }
        }
        stateChecks++;
    }
    private static void queries(Sketch a, Sketch b, List<Double> xs) throws Exception {
        List<Integer> x = a.query(xs), y = b.query(xs);
        if (!x.equals(y)) throw new AssertionError("rank mismatch");
        rankChecks += xs.size(); state(a, b);
    }
    private static boolean transition(Sketch a, Sketch b, Runnable first, Runnable second) throws Exception {
        Throwable x = null, y = null;
        try { first.run(); } catch (RuntimeException | AssertionError error) { x = error; }
        try { second.run(); } catch (RuntimeException | AssertionError error) { y = error; }
        if (x == null && y == null) return true;
        if (x == null || y == null || x.getClass() != y.getClass() || !Objects.equals(x.getMessage(),y.getMessage()))
            throw new AssertionError("mutation failure mismatch",x == null ? y : x);
        state(a,b); sharedFailures++;
        System.err.println("Shared upstream mutation failure: " + x);
        return false;
    }
    private static void kernels() {
        for (int n : new int[]{3,4,8,32,128,1024}) for (int seed = 0; seed < 100; seed++) {
            RandomBits random = new RandomBits(seed);
            double[] x = new double[n]; int[] y = new int[n];
            for (int i = 0; i < n; i++) {
                x[i] = (i == 0 ? -100 : x[i - 1]) + (1 + random.bounded(1000000)) / 1000000.0;
                y[i] = i == 0 ? 0 : y[i - 1] + (int)random.bounded(100000);
            }
            kernel(x, y);
        }
        kernel(new double[]{0, Double.MIN_VALUE, 2 * Double.MIN_VALUE, 3 * Double.MIN_VALUE}, new int[]{0,1,2,3});
        kernel(new double[]{-Double.MAX_VALUE,-1,0,1,Double.MAX_VALUE}, new int[]{0,1,2,3,4});
        kernel(new double[]{-2,-1,-0.0,0.0,1,2}, new int[]{0,1,1,1,1,2});
        kernel(new double[]{-2,-1,0,1,2}, new int[]{Integer.MIN_VALUE,Integer.MAX_VALUE,0,Integer.MAX_VALUE,Integer.MIN_VALUE});
    }
    private static void kernel(double[] x, int[] y) {
        double[][] derivatives = {
            new OriginalSplineSketch.PchipInterpolator(x,y).getDerivatives(), new SplineSketch.PchipInterpolator(x,y).getDerivatives(),
            new OriginalSplineSketchMG.PchipInterpolator(x,y).getDerivatives(), new SplineSketchMG.PchipInterpolator(x,y).getDerivatives()
        };
        for (int j : new int[]{0,2}) for (int i = 0; i < x.length; i++) {
            if (Double.doubleToLongBits(derivatives[j][i]) != Double.doubleToLongBits(derivatives[j+1][i]))
                throw new AssertionError("derivative mismatch");
            derivativeChecks++;
        }
    }
    private static void checks() throws Exception {
        kernels();
        // Retain both signed zeros in MG and force its indexed large-batch path.
        Sketch zeroA = new Sketch(true,true,32), zeroB = new Sketch(false,true,32);
        double[] zeros = new double[200];
        for (int i = 0; i < zeros.length; i++) {
            zeros[i] = i % 20 == 0 ? -0.0 : i % 20 == 1 ? 0.0 : i % 20 - 1;
            zeroA.update(zeros[i]); zeroB.update(zeros[i]);
        }
        zeroA.finish(); zeroB.finish(); queries(zeroA,zeroB,probes(zeros));
        int cases = 0;
        for (boolean mg : new boolean[]{false,true}) for (int k : new int[]{6,8,32,128,512})
        for (int shape = 0; shape < 5; shape++) for (int seed : new int[]{0,1,2,71,72,73}) {
            double[] input = data(shape, seed, 6000); List<Double> xs = probes(input);
            Sketch a = new Sketch(true,mg,k), b = new Sketch(false,mg,k);
            queries(a,b,xs);
            for (int i = 0; i < input.length; i++) {
                a.update(input[i]); b.update(input[i]);
                if ((i + 1) % (5 * k) == 0 || i == 0 || i == 5 * k - 2 || i == 5 * k || i == input.length - 2) {
                    queries(a,b,xs);
                    queries(a,b,Collections.singletonList(input[i]));
                    queries(a,b,Collections.emptyList());
                }
            }
            a.finish(); b.finish(); queries(a,b,xs);
            // Reverse/repeated queries exercise the MG batch index independently
            // of query order, including both signed zeros and infinities.
            List<Double> reversed = new ArrayList<>(xs); Collections.reverse(reversed); queries(a,b,reversed);
            boolean completed = true;
            for (int newK : new int[]{Math.max(6,k/2),2*k,k}) {
                if (!transition(a,b,() -> a.resize(newK),() -> b.resize(newK))) { completed = false; break; }
                queries(a,b,xs);
            }
            if (completed && mg) {
                if (transition(a,b,a::compress,b::compress)) queries(a,b,xs);
            } else if (completed) {
                Sketch c = new Sketch(true,false,k), d = new Sketch(false,false,k);
                for (double v : data(shape,seed + 17,6000)) { c.update(v); d.update(v); }
                c.finish(); d.finish();
                Object u = OriginalSplineSketch.merge((OriginalSplineSketch)a.value,(OriginalSplineSketch)c.value);
                Object v = SplineSketch.merge((SplineSketch)b.value,(SplineSketch)d.value);
                same(((OriginalSplineSketch)u).query(xs), ((SplineSketch)v).query(xs));
                rankChecks += xs.size(); state(a,b); state(c,d); mergeChecks++;
            }
            cases++;
        }
        System.out.printf(Locale.ROOT,"{\"cases\":%d,\"rank_queries\":%d,\"full_state_checks\":%d,\"derivative_values\":%d,\"plain_merges\":%d,\"shared_upstream_failures\":%d}%n",
                          cases,rankChecks,stateChecks,derivativeChecks,mergeChecks,sharedFailures);
    }

    private static final Map<Class<?>, List<Field>> graphFields = new HashMap<>();
    private static long bytes(Object root) throws Exception {
        Set<Object> seen = Collections.newSetFromMap(new IdentityHashMap<>());
        ArrayDeque<Object> todo = new ArrayDeque<>(); todo.add(root); long total = 0;
        while (!todo.isEmpty()) {
            Object object = todo.removeLast(); if (!seen.add(object)) continue;
            total += instrumentation.getObjectSize(object);
            if (object.getClass().isArray()) {
                if (!object.getClass().getComponentType().isPrimitive())
                    for (Object item : (Object[])object) if (item != null) todo.add(item);
            } else {
                List<Field> fields = graphFields.get(object.getClass());
                if (fields == null) {
                    fields = new ArrayList<>();
                    for (Class<?> type = object.getClass(); type != null; type = type.getSuperclass())
                        for (Field field : type.getDeclaredFields())
                            if (!Modifier.isStatic(field.getModifiers()) && !field.getType().isPrimitive()) {
                                field.setAccessible(true); fields.add(field);
                            }
                    graphFields.put(object.getClass(),fields);
                }
                for (Field field : fields) { Object item = field.get(object); if (item != null) todo.add(item); }
            }
        }
        return total;
    }
    private static long[] timeQuery(Sketch sketch, List<Double> xs, int passes, com.sun.management.ThreadMXBean bean) {
        long thread = Thread.currentThread().threadId(), allocation = bean.getThreadAllocatedBytes(thread), start = System.nanoTime();
        long sum = 0;
        for (int pass = 0; pass < passes; pass++) for (int rank : sketch.query(xs)) sum += rank;
        long elapsed = System.nanoTime() - start, allocated = bean.getThreadAllocatedBytes(thread) - allocation;
        sink += sum; return new long[]{elapsed, allocated};
    }
    private static void benchmark(int trials) throws Exception {
        if (instrumentation == null) throw new IllegalArgumentException("benchmark requires the agent");
        var bean = (com.sun.management.ThreadMXBean)ManagementFactory.getThreadMXBean();
        if (!bean.isThreadAllocatedMemorySupported()) throw new IllegalStateException("thread allocation measurement unavailable");
        bean.setThreadAllocatedMemoryEnabled(true);
        System.out.println("variant,mg,k,shape,seed,trial,resident_bytes,update_ns,update_allocated_bytes,batch_ns,batch_allocated_bytes,single_ns,single_allocated_bytes,query_hash");
        for (boolean mg : new boolean[]{false,true}) for (int k : new int[]{8,32,128,512,1024}) {
            double[] warm = data(3,70,6000); List<Double> warmQueries = probes(warm);
            for (boolean original : new boolean[]{true,false}) for (int pass = 0; pass < 4; pass++) {
                Sketch sketch = new Sketch(original,mg,k); for (double x : warm) sketch.update(x); sketch.finish();
                timeQuery(sketch,warmQueries,8,bean);
                for (int i = 0; i < 128; i++) timeQuery(sketch,Collections.singletonList(warmQueries.get(i)),4,bean);
            }
            for (int shape = 0; shape < 5; shape++) for (int seed : new int[]{71,72,73}) {
                double[] input = data(shape,seed,6000); List<Double> xs = probes(input);
                List<List<Double>> singletons = new ArrayList<>();
                for (int i = 0; i < 128; i++) singletons.add(Collections.singletonList(xs.get(i * (xs.size()-1)/127)));
                for (int trial = 0; trial < trials; trial++) for (int position = 0; position < 2; position++) {
                    boolean original = ((trial + position) & 1) == 0;
                    Sketch sketch = new Sketch(original,mg,k);
                    long thread = Thread.currentThread().threadId(), allocated = bean.getThreadAllocatedBytes(thread), start = System.nanoTime();
                    for (double x : input) sketch.update(x); sketch.finish();
                    long update = System.nanoTime() - start, updateAllocated = bean.getThreadAllocatedBytes(thread) - allocated;
                    long retained = bytes(sketch.value);
                    List<Integer> expected = sketch.query(xs); long hash = 0xcbf29ce484222325L;
                    for (int rank : expected) hash = (hash ^ rank) * 0x100000001b3L;
                    timeQuery(sketch,xs,4,bean); long[] batch = timeQuery(sketch,xs,16,bean);
                    long singleTime = 0, singleAllocated = 0;
                    for (List<Double> one : singletons) {
                        long[] result = timeQuery(sketch,one,8,bean); singleTime += result[0]; singleAllocated += result[1];
                    }
                    System.out.printf(Locale.ROOT,"%s,%s,%d,%d,%d,%d,%d,%.9f,%d,%.9f,%.9f,%.9f,%.9f,%s%n",
                        original?"original":"optimized",mg,k,shape,seed,trial,retained,update/6000.0,updateAllocated,
                        batch[0]/(16.0*xs.size()),batch[1]/(16.0*xs.size()),singleTime/1024.0,singleAllocated/1024.0,Long.toUnsignedString(hash));
                }
            }
            System.out.flush();
            System.err.println("Benchmark complete: mg=" + mg + " k=" + k);
        }
    }
    public static void main(String[] args) throws Exception {
        if (args.length == 1 && args[0].equals("--check")) checks();
        else if (args.length == 2 && args[0].equals("--benchmark")) benchmark(Integer.parseInt(args[1]));
        else throw new IllegalArgumentException("usage: --check | --benchmark TRIALS");
    }
}
