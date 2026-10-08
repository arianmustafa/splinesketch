import java.lang.instrument.Instrumentation;
import java.lang.reflect.*;
import java.nio.file.*;
import java.util.*;

/** Adapter for the unmodified pinned authors' classes. No third-party jars. */
public class AuthorsComparison {
    private static Instrumentation instrumentation;
    public static void premain(String options, Instrumentation inst) { instrumentation = inst; }
    private static final Map<Class<?>, List<Field>> fields = new HashMap<>();
    private static List<Field> references(Class<?> cls) {
        return fields.computeIfAbsent(cls, key -> {
            List<Field> result = new ArrayList<>();
            for (Class<?> c=key; c!=null; c=c.getSuperclass()) for (Field f:c.getDeclaredFields()) {
                if (!Modifier.isStatic(f.getModifiers()) && !f.getType().isPrimitive()) {
                    f.setAccessible(true); result.add(f);
                }
            }
            return result;
        });
    }
    // Object graph rooted at the sketch, with identity deduplication. Includes
    // arrays, map table/nodes, boxed keys and counters, and object headers.
    private static long bytes(Object root) throws Exception {
        Set<Object> seen = Collections.newSetFromMap(new IdentityHashMap<>());
        ArrayDeque<Object> queue = new ArrayDeque<>(); queue.add(root);
        long total=0;
        while (!queue.isEmpty()) {
            Object o=queue.removeLast(); if(!seen.add(o)) continue;
            total += instrumentation.getObjectSize(o);
            if (o.getClass().isArray()) {
                if (!o.getClass().getComponentType().isPrimitive())
                    for(Object x:(Object[])o) if(x!=null) queue.add(x);
            } else for(Field f:references(o.getClass())) { Object x=f.get(o); if(x!=null) queue.add(x); }
        }
        return total;
    }
    private static List<Double> read(Path path) throws Exception {
        List<Double> result=new ArrayList<>();
        for(String line:Files.readAllLines(path)) result.add(Double.parseDouble(line));
        return result;
    }
    private static class Sketch {
        final Object value; final boolean mg;
        Sketch(int k,boolean mg) { this.mg=mg; value=mg ? new SplineSketchMG(k) : new SplineSketch(k); }
        void update(double x) { if(mg) ((SplineSketchMG)value).update(x); else ((SplineSketch)value).update(x); }
        void finish() { if(mg) ((SplineSketchMG)value).consolidate(); else ((SplineSketch)value).consolidate(); }
        List<Integer> query(List<Double> x) { return mg ? ((SplineSketchMG)value).query(x) : ((SplineSketch)value).query(x); }
        int compactEstimate() { return mg ? ((SplineSketchMG)value).serializedSketchBytesCompact() : ((SplineSketch)value).serializedSketchBytesCompact(); }
    }
    private static int upperBound(List<Double> sorted,double x) {
        int lo=0,hi=sorted.size(); while(lo<hi) {int mid=(lo+hi)>>>1;if(sorted.get(mid)<=x)lo=mid+1;else hi=mid;}return lo;
    }
    private static volatile long sink;
    public static void main(String[] args) throws Exception {
        if((args.length!=1 && args.length!=2) || instrumentation==null) throw new IllegalArgumentException("requires agent and data directory");
        System.out.println("variant,k,shape,seed,resident_peak,final_bytes,compact_estimate,update_allocated_bytes,update_ns,rank_batch_ns,rank_single_ns,median_percent,p95_percent,max_percent");
        var bean=(com.sun.management.ThreadMXBean)java.lang.management.ManagementFactory.getThreadMXBean();
        long thread=Thread.currentThread().threadId();
        int[] capacities=args.length==2 ? Arrays.stream(args[1].split(",")).mapToInt(Integer::parseInt).toArray() : new int[]{6,8,12,16,24,32,48,64,96,128,192,256};
        for(boolean mg:new boolean[]{false,true}) for(int k:capacities) {
            // Warm the same code path before timed trials; graph measurement is
            // a separate replay and never included in timing or allocation totals.
            List<Double> warm=read(Path.of(args[0],"outliers-0.data"));
            for(int r=0;r<5;r++){Sketch s=new Sketch(k,mg);for(double x:warm)s.update(x);s.finish();sink+=s.query(warm).get(0);}
            for(String shape:new String[]{"uniform","clustered","duplicates","outliers","ascending"}) for(int seed=0;seed<3;seed++) {
                String name=shape+"-"+seed;
                List<Double> data=read(Path.of(args[0],name+".data"));
                List<Double> queries=read(Path.of(args[0],name+".queries"));
                Sketch measured=new Sketch(k,mg);long peak=bytes(measured.value);
                // Retained fields only change at consolidation. update() merely
                // writes into the preallocated primitive buffer otherwise.
                Field bufferIndex=measured.value.getClass().getDeclaredField("bufferIndex"); bufferIndex.setAccessible(true);
                for(double x:data) {measured.update(x);if(bufferIndex.getInt(measured.value)==0)peak=Math.max(peak,bytes(measured.value));}
                measured.finish();long end=bytes(measured.value);peak=Math.max(peak,end);
                Sketch s=new Sketch(k,mg);
                long alloc0=bean.getThreadAllocatedBytes(thread),t0=System.nanoTime();
                for(double x:data)s.update(x);s.finish();
                long t1=System.nanoTime(),allocated=bean.getThreadAllocatedBytes(thread)-alloc0;
                List<Integer> ranks=s.query(queries);long t2=System.nanoTime();
                for(double x:queries)sink+=s.query(Collections.singletonList(x)).get(0);
                long t3=System.nanoTime();
                List<Double> sorted=new ArrayList<>(data);Collections.sort(sorted);
                double[] errors=new double[queries.size()];
                for(int i=0;i<errors.length;i++) {
                    int rank=ranks.get(i);if(rank<0||rank>data.size())throw new AssertionError("invalid rank");
                    errors[i]=100.0*Math.abs(rank-upperBound(sorted,queries.get(i)))/data.size();
                }
                Arrays.sort(errors);
                System.out.printf(Locale.ROOT,"%s,%d,%s,%d,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.9f,%.9f,%.9f%n",
                    mg?"java_mg":"java",k,shape,seed,peak,end,s.compactEstimate(),allocated,
                    (t1-t0)/(double)data.size(),(t2-t1)/(double)queries.size(),(t3-t2)/(double)queries.size(),
                    errors[(errors.length-1)/2],errors[(errors.length-1)*95/100],errors[errors.length-1]);
            }
        }
    }
}
